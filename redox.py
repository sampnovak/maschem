"""Oxidation-state assignment and redox bookkeeping."""

from __future__ import annotations

from dataclasses import dataclass, field

from .data import E0
from .elements import ALKALI, ALKALINE, ELEMENTS, HALOGENS, METALS
from .formula import Species, hill_formula

# Species whose oxidation states cannot be derived by the simple rules.
_OVERRIDES: dict[str, dict[str, float]] = {
    "H2O2": {"H": 1, "O": -1},
    "Na2O2": {"Na": 1, "O": -1},
    "BaO2": {"Ba": 2, "O": -1},
    "CaO2": {"Ca": 2, "O": -1},
    "KO2": {"K": 1, "O": -0.5},
    "O3": {"O": 0},
    "OF2": {"O": 2, "F": -1},
    "FeS2": {"Fe": 2, "S": -1},
    "Fe3O4": {"Fe": 8 / 3, "O": -2},
    "Mn3O4": {"Mn": 8 / 3, "O": -2},
    "Pb3O4": {"Pb": 8 / 3, "O": -2},
    "N2H4": {"N": -2, "H": 1},
    "NH4NO3": {"N": 1, "H": 1, "O": -2},
    "S8": {"S": 0},
    "P4": {"P": 0},
    "C6H12O6": {"C": 0, "H": 1, "O": -2},
}

_DEFAULT_NEGATIVE = {
    "O": -2, "S": -2, "Se": -2, "Te": -2,
    "N": -3, "P": -3, "As": -3,
    "C": -4, "Si": -4, "B": -3, "H": -1,
}


def potential(couple: str) -> float | None:
    """Standard reduction potential for an ``Ox/Red`` couple, volts vs SHE."""
    return E0.get(couple)


def cell_potential(cathode: str, anode: str) -> float | None:
    """E°cell when ``cathode`` oxidises ``anode``; positive means favourable."""
    e_cathode = E0.get(cathode)
    e_anode = E0.get(anode)
    if e_cathode is None or e_anode is None:
        return None
    return e_cathode - e_anode


def oxidation_states(sp: Species) -> dict[str, float]:
    """Assign an oxidation state to each element in a species."""
    counts = sp.counts
    if len(counts) == 1:
        symbol = next(iter(counts))
        return {symbol: sp.charge / counts[symbol]}
    key = hill_formula(counts)
    if key in _OVERRIDES:
        states = dict(_OVERRIDES[key])
        # scale nothing; overrides are per atom
        return states

    states: dict[str, float] = {}

    if "F" in counts:
        states["F"] = -1.0
        if "O" in counts:
            states["O"] = 2.0
    elif "O" in counts:
        oxygen = -2.0
        if counts.get("O", 0) == 2 and set(counts) <= {"H", "O"} and counts.get("H", 0) == 2:
            oxygen = -1.0
        states["O"] = oxygen

    for symbol in counts:
        if symbol in ALKALI:
            states[symbol] = 1.0
        elif symbol in ALKALINE:
            states[symbol] = 2.0
    if "Al" in counts:
        states["Al"] = 3.0
    if "Zn" in counts:
        states["Zn"] = 2.0
    if "Ag" in counts:
        states["Ag"] = 1.0

    if "H" in counts and "H" not in states:
        others = [s for s in counts if s != "H"]
        metallic = others and all(s in METALS or s == "B" for s in others)
        states["H"] = -1.0 if metallic else 1.0

    if "O" not in counts and "F" not in counts:
        for symbol in counts:
            if symbol in HALOGENS and symbol not in states:
                states[symbol] = -1.0

    def residual() -> float:
        used = sum(states[s] * counts[s] for s in states)
        return sp.charge - used

    remaining = [s for s in counts if s not in states]

    if len(remaining) == 1:
        symbol = remaining[0]
        states[symbol] = residual() / counts[symbol]
        return states

    # Give metals their most common positive state first.
    for symbol in list(remaining):
        element = ELEMENTS.get(symbol)
        if symbol in METALS and element and element.states:
            positive = next((s for s in element.states if s > 0), None)
            if positive is not None:
                states[symbol] = float(positive)
                remaining.remove(symbol)

    if not remaining:
        return states

    # All but the most electronegative leftover get a typical negative state.
    remaining.sort(key=lambda s: ELEMENTS[s].electronegativity if s in ELEMENTS else 0.0)
    for symbol in remaining[:-1]:
        states[symbol] = float(_DEFAULT_NEGATIVE.get(symbol, -1))

    last = remaining[-1]
    states[last] = residual() / counts[last]
    return states


@dataclass
class RedoxChange:
    element: str
    before: float
    after: float
    atoms: float
    electrons: float
    oxidised: bool
    from_species: str = ""
    to_species: str = ""

    @property
    def per_atom(self) -> float:
        return self.after - self.before

    @property
    def transition(self) -> str:
        """Species on each side, when the element is carried by two of them."""
        if not self.from_species and not self.to_species:
            return ""
        return f"{self.from_species or '?'} \u2192 {self.to_species or '?'}"


@dataclass
class RedoxResult:
    happened: bool
    electrons: float = 0.0
    changes: list[RedoxChange] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def oxidised_elements(self) -> list[str]:
        return list(dict.fromkeys(c.element for c in self.changes if c.oxidised))

    def reduced_elements(self) -> list[str]:
        return list(dict.fromkeys(c.element for c in self.changes if not c.oxidised))


def analyse(
    reactants: list[Species],
    rcoeffs: list[int],
    products: list[Species],
    pcoeffs: list[int],
) -> RedoxResult:
    """Find which atoms change oxidation state, one species at a time.

    States are collected per species rather than per element, so an element
    carried by several species keeps them apart: in KMnO4 + HCl the chlorine
    of KCl/MnCl2 is still -1 while the chlorine of Cl2 is 0, instead of both
    being averaged into a meaningless -0.38.
    """
    left: dict[str, list[list]] = {}
    right: dict[str, list[list]] = {}
    for side, species_list, coeffs in (
        (left, reactants, rcoeffs),
        (right, products, pcoeffs),
    ):
        for sp, coeff in zip(species_list, coeffs):
            if coeff == 0:
                continue
            for symbol, state in oxidation_states(sp).items():
                atoms = sp.counts.get(symbol, 0.0) * coeff
                if atoms > 0:
                    side.setdefault(symbol, []).append([state, atoms, sp.formula])

    changes: list[RedoxChange] = []
    for symbol in sorted(set(left) | set(right)):
        left_entries = [entry[:] for entry in left.get(symbol, [])]
        right_entries = [entry[:] for entry in right.get(symbol, [])]
        if not left_entries or not right_entries:
            continue
        # Atoms that keep their oxidation state are not part of the change.
        for before in left_entries:
            for after in right_entries:
                if abs(before[0] - after[0]) < 1e-9 and before[1] > 0 and after[1] > 0:
                    shared = min(before[1], after[1])
                    before[1] -= shared
                    after[1] -= shared
        left_entries = [e for e in left_entries if e[1] > 1e-12]
        right_entries = [e for e in right_entries if e[1] > 1e-12]
        while left_entries and right_entries:
            before = max(left_entries, key=lambda e: e[1])
            after = max(right_entries, key=lambda e: e[1])
            atoms = min(before[1], after[1])
            delta = (after[0] - before[0]) * atoms
            changes.append(
                RedoxChange(
                    element=symbol,
                    before=before[0],
                    after=after[0],
                    atoms=atoms,
                    electrons=delta,
                    oxidised=delta > 0,
                    from_species=before[2],
                    to_species=after[2],
                )
            )
            before[1] -= atoms
            after[1] -= atoms
            if before[1] <= 1e-12:
                left_entries.remove(before)
            if after[1] <= 1e-12:
                right_entries.remove(after)

    if not changes:
        return RedoxResult(happened=False)

    changes.sort(key=lambda c: -abs(c.electrons))
    electrons = sum(c.electrons for c in changes if c.electrons > 0)
    return RedoxResult(happened=True, electrons=electrons, changes=changes)


def oxidation_table(species: list[Species]) -> list[tuple[str, dict[str, float]]]:
    return [(sp.display, oxidation_states(sp)) for sp in species]
