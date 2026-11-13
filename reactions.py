"""Product prediction: known reactions first, then inorganic heuristics."""

from __future__ import annotations

import re
from dataclasses import dataclass

from .elements import (
    ALKALI,
    ELEMENTS,
    HALOGENS,
    METALS,
    NONMETALS,
    REACTIVITY_SERIES,
)
from .formula import (
    FormulaError,
    Species,
    formula_from_ions,
    match_key,
    parse_formula,
    parse_species,
)
from .redox import cell_potential

ANION_CHARGES: dict[str, int] = {
    "F": -1, "Cl": -1, "Br": -1, "I": -1,
    "NO3": -1, "NO2": -1, "ClO3": -1, "ClO4": -1, "ClO": -1,
    "SO4": -2, "SO3": -2, "CO3": -2, "SiO3": -2, "CrO4": -2,
    "Cr2O7": -2, "C2O4": -2, "HPO4": -2,
    "PO4": -3, "PO3": -3, "N": -3, "P": -3,
    "OH": -1, "CH3COO": -1, "C2H3O2": -1, "HCOO": -1, "CN": -1, "SCN": -1,
    "HCO3": -1, "HSO4": -1, "MnO4": -1, "S": -2, "O": -2, "Se": -2,
}

CATION_CHARGES: dict[str, int] = {
    "NH4": 1,
    "Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1,
    "Be": 2, "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2,
    "Al": 3, "Zn": 2, "Ag": 1, "Cu": 2, "Fe": 2, "Pb": 2, "Sn": 2,
    "Ni": 2, "Co": 2, "Mn": 2, "Cr": 3, "Hg": 2, "Cd": 2, "Bi": 3,
}

_SALT_CACHE: dict[str, tuple[str, int, str, int] | None] = {}


def salt_parts(sp: Species) -> tuple[str, int, str, int] | None:
    """Decompose a neutral salt into (cation, charge, anion, charge)."""
    key = sp.formula
    if key in _SALT_CACHE:
        return _SALT_CACHE[key]
    result = _split_salt_uncached(sp)
    _SALT_CACHE[key] = result
    return result


def _split_salt_uncached(sp: Species) -> tuple[str, int, str, int] | None:
    counts = {sym: n for sym, n in sp.counts.items()}
    if sp.charge != 0:
        return None

    cation_candidates: list[tuple[str, int, dict[str, float]]] = []
    for symbol in sorted(METALS):
        if symbol in counts:
            # Any of the common positive states can be the cation here, so
            # FeCl3 and FeCl2 are both recognised.
            for state in ELEMENTS[symbol].states:
                if state > 0:
                    cation_candidates.append((symbol, state, {symbol: 1.0}))
    if "N" in counts and "H" in counts and counts.get("H", 0) >= 4:
        cation_candidates.append(("NH4", 1, {"N": 1.0, "H": 4.0}))

    for cat, cat_charge, cat_comp in cation_candidates:
        for anion, an_charge in ANION_CHARGES.items():
            if anion == "O":
                continue  # a metal oxide is not a salt, even though it splits this way
            an_comp = _composition_of(anion)
            if not an_comp:
                continue
            for n_cat in range(1, 9):
                for n_an in range(1, 9):
                    combined: dict[str, float] = {}
                    for sym, n in cat_comp.items():
                        combined[sym] = combined.get(sym, 0) + n * n_cat
                    for sym, n in an_comp.items():
                        combined[sym] = combined.get(sym, 0) + n * n_an
                    if combined != counts:
                        continue
                    if n_cat * cat_charge + n_an * an_charge == sp.charge:
                        return (cat, cat_charge, anion, an_charge)
    return None


_COMPOSITION_CACHE: dict[str, dict[str, float] | None] = {}


def _composition_of(formula: str) -> dict[str, float] | None:
    if formula in _COMPOSITION_CACHE:
        return _COMPOSITION_CACHE[formula]
    try:
        comp = parse_formula(formula)
    except FormulaError:
        comp = None
    _COMPOSITION_CACHE[formula] = comp
    return comp


def is_soluble(cation: str, anion: str) -> bool:
    from .data import (
        INSOLUBLE_ANIONS,
        INSOLUBLE_HALIDES,
        INSOLUBLE_SULFATES,
        SOLUBLE_HYDROXIDES,
        SOLUBLE_ALWAYS,
    )

    if anion in SOLUBLE_ALWAYS:
        return True
    if anion in ("Cl", "Br", "I"):
        return cation not in INSOLUBLE_HALIDES
    if anion == "SO4":
        return cation not in INSOLUBLE_SULFATES
    if anion == "OH":
        return cation in SOLUBLE_HYDROXIDES
    if anion in INSOLUBLE_ANIONS:
        return cation == "NH4" or cation in ALKALI
    return True


# ---------------------------------------------------------------------------
class PredictionError(ValueError):
    pass


@dataclass
class Prediction:
    products: list[Species]
    source: str
    notes: list[str]


@dataclass(frozen=True)
class _Entry:
    reactants: tuple[str, ...]
    products: tuple[str, ...]
    note: str


_ENTRIES: list[_Entry] = []


def _add(reactants, products, note):
    _ENTRIES.append(_Entry(tuple(reactants), tuple(products), note))


_add(("H2", "O2"), ("H2O(l)",), "hydrogen burns in oxygen; very exothermic")
_add(("CH4", "O2"), ("CO2", "H2O"), "complete combustion of methane")
_add(("C2H6", "O2"), ("CO2", "H2O"), "complete combustion of ethane")
_add(("C3H8", "O2"), ("CO2", "H2O"), "complete combustion of propane")
_add(("C4H10", "O2"), ("CO2", "H2O"), "complete combustion of butane")
_add(("C2H5OH", "O2"), ("CO2", "H2O"), "combustion of ethanol")
_add(("C6H6", "O2"), ("CO2", "H2O"), "combustion of benzene (sooty flame)")
_add(("C", "O2"), ("CO2",), "carbon burns in excess oxygen")
_add(("S", "O2"), ("SO2",), "sulfur burns with a blue flame")
_add(("N2", "O2"), ("NO",), "needs very high temperature (lightning, engines)")
_add(("Fe", "S"), ("FeS",), "iron and sulfur glow when heated together")
_add(("Zn", "S"), ("ZnS",), "zinc blende formation")
_add(("Cu", "S"), ("CuS",), "copper(II) sulfide")
_add(("Al", "S"), ("Al2S3",), "aluminium sulfide")
_add(("Na", "Cl2"), ("NaCl",), "sodium burns in chlorine")
_add(("H2", "Cl2"), ("HCl",), "photochemical combination")
_add(("H2", "N2"), ("NH3",), "Haber process: catalyst, high pressure, ~700 K")
_add(("Mg", "O2"), ("MgO",), "brilliant white flame")
_add(("Al", "O2"), ("Al2O3",), "protective oxide layer")
_add(("Fe", "O2"), ("Fe2O3",), "rusting is slow; also needs water")
_add(("Cu", "O2"), ("CuO",), "copper tarnishes when heated")
_add(("Zn", "O2"), ("ZnO",), "zinc oxide")
_add(("Na", "O2"), ("Na2O2",), "sodium peroxide (excess oxygen)")
_add(("Ca", "O2"), ("CaO",), "calcium oxide")
_add(("CaO", "H2O"), ("Ca(OH)2",), "slaking of lime; strongly exothermic")
_add(("Na2O", "H2O"), ("NaOH",), "metal oxide plus water gives hydroxide")
_add(("K2O", "H2O"), ("KOH",), "metal oxide plus water gives hydroxide")
_add(("SO3", "H2O"), ("H2SO4",), "sulfur trioxide plus water gives sulfuric acid")
_add(("SO2", "H2O"), ("H2SO3",), "sulfurous acid (weak)")
_add(("CO2", "H2O"), ("H2CO3",), "carbonic acid (weak, equilibrium)")
_add(("N2O5", "H2O"), ("HNO3",), "nitric acid from its anhydride")
_add(("P4O10", "H2O"), ("H3PO4",), "phosphoric acid from phosphorus pentoxide")
_add(("Na", "H2O"), ("NaOH", "H2"), "violent; hydrogen ignites")
_add(("K", "H2O"), ("KOH", "H2"), "violent; lilac flame, hydrogen ignites")
_add(("Li", "H2O"), ("LiOH", "H2"), "steady reaction, hydrogen bubbles")
_add(("Ca", "H2O"), ("Ca(OH)2", "H2"), "calcium granules fizz in water")
_add(("Mg", "H2O"), ("Mg(OH)2", "H2"), "slow in cold water, fast in steam")
_add(("Fe", "HCl"), ("FeCl2", "H2"), "dilute acid; gives iron(II)")
_add(("Zn", "HCl"), ("ZnCl2", "H2"), "bubbles of hydrogen")
_add(("Mg", "HCl"), ("MgCl2", "H2"), "rapid hydrogen evolution")
_add(("Al", "HCl"), ("AlCl3", "H2"), "needs the oxide layer removed")
_add(("NaOH", "HCl"), ("NaCl", "H2O"), "neutralisation, \u0394H \u2248 -57 kJ/mol")
_add(("KOH", "H2SO4"), ("K2SO4", "H2O"), "neutralisation")
_add(("Ca(OH)2", "HCl"), ("CaCl2", "H2O"), "neutralisation")
_add(("NH3", "HCl"), ("NH4Cl",), "white smoke of ammonium chloride")
_add(("NaHCO3", "HCl"), ("NaCl", "H2O", "CO2"), "baking soda plus acid: fizzing")
_add(("CaCO3", "HCl"), ("CaCl2", "H2O", "CO2"), "limestone dissolves with fizzing")
_add(("Na2CO3", "HCl"), ("NaCl", "H2O", "CO2"), "carbonate plus acid")
_add(("CaCO3",), ("CaO", "CO2"), "thermal decomposition (calcination, ~1200 K)")
_add(("KClO3",), ("KCl", "O2"), "needs heat and a MnO2 catalyst")
_add(("H2O2",), ("H2O", "O2"), "slowly decomposes; fast with MnO2 or catalase")
_add(("Fe2O3", "Al"), ("Al2O3", "Fe"), "thermite: reaches ~2500 \u00b0C, molten iron")
_add(("Fe", "CuSO4"), ("FeSO4", "Cu"), "iron displaces copper; blue colour fades")
_add(("Zn", "CuSO4"), ("ZnSO4", "Cu"), "zinc displaces copper")
_add(("Cu", "AgNO3"), ("Cu(NO3)2", "Ag"), "copper displaces silver; silver crystals form")
_add(("AgNO3", "NaCl"), ("AgCl", "NaNO3"), "white precipitate of silver chloride")
_add(("BaCl2", "Na2SO4"), ("BaSO4", "NaCl"), "white precipitate of barium sulfate")
_add(("Pb(NO3)2", "KI"), ("PbI2", "KNO3"), "golden-yellow lead iodide precipitate")
_add(("CaCl2", "Na2CO3"), ("CaCO3", "NaCl"), "white calcium carbonate precipitate")
_add(("CuSO4", "NaOH"), ("Cu(OH)2", "Na2SO4"), "pale blue copper hydroxide precipitate")
_add(("FeCl3", "NaOH"), ("Fe(OH)3", "NaCl"), "red-brown iron(III) hydroxide precipitate")
_add(("AgNO3", "NaBr"), ("AgBr", "NaNO3"), "cream silver bromide precipitate")
_add(("BaCl2", "H2SO4"), ("BaSO4", "HCl"), "white barium sulfate precipitate")
_add(("NH4Cl", "NaOH"), ("NH3", "H2O", "NaCl"), "releases ammonia gas when warmed")
_add(("Mg", "CO2"), ("MgO", "C"), "burning magnesium keeps burning in CO2")
_add(("CuO", "H2"), ("Cu", "H2O"), "hydrogen reduces hot copper(II) oxide")
_add(("Fe2O3", "CO"), ("Fe", "CO2"), "blast-furnace reduction")
_add(("CuO", "C"), ("Cu", "CO2"), "carbon reduces copper(II) oxide")
_add(("Cl2", "NaBr"), ("NaCl", "Br2"), "chlorine displaces bromine")
_add(("Cl2", "KI"), ("KCl", "I2"), "chlorine displaces iodine")
_add(("CH3COOH", "NaHCO3"), ("CH3COONa", "H2O", "CO2"), "vinegar plus baking soda")


def _entry_keys(entry: _Entry) -> set[str]:
    keys = set()
    for text in entry.reactants:
        try:
            keys.add(match_key(parse_species(text)))
        except FormulaError:
            keys.add(text)
    return keys


_ENTRY_TABLE = [(_entry_keys(e), e) for e in _ENTRIES]


def _is_element(sp: Species) -> bool:
    return len(sp.counts) == 1


def _element_symbol(sp: Species) -> str | None:
    if _is_element(sp):
        return next(iter(sp.counts))
    return None


def _anion_charge(symbol: str) -> int | None:
    if symbol in HALOGENS:
        return -1
    return {
        "O": -2, "S": -2, "Se": -2, "Te": -2,
        "N": -3, "P": -3, "As": -3,
        "C": -4, "Si": -4, "B": -3, "H": -1,
    }.get(symbol)


def _acid_parts(sp: Species) -> tuple[str, int, int] | None:
    """Return (anion formula, anion charge, acidic protons) for a neutral acid."""
    if sp.charge != 0:
        return None
    counts = sp.counts
    hydrogens = int(counts.get("H", 0))
    if hydrogens <= 0 or len(counts) < 2:
        return None
    if set(counts) == {"H", "O"}:
        return None  # water and hydrogen peroxide are not acids here
    return _anion_of_acid(sp.formula, hydrogens)


_LEADING_H = re.compile(r"^H(\d*)")


def _anion_of_acid(formula: str, hydrogens: int) -> tuple[str, int, int] | None:
    """(anion formula, anion charge, acidic protons) for an acid.

    Works on the typed formula so the result keeps its conventional order
    (``H2SO4`` gives ``SO4``, ``CH3COOH`` gives ``CH3COO``) instead of a
    Hill-sorted one.  Only recognised anions qualify, which keeps hydroxides
    and other bases from being mistaken for acids.
    """
    candidates: list[tuple[str, int]] = []
    match = _LEADING_H.match(formula)
    if match and int(match.group(1) or 1) == hydrogens:
        body = formula[match.end() :]
        if body:
            candidates.append((body, hydrogens))
    if formula.endswith("COOH"):
        # A carboxylic acid gives up one proton; the other H atoms stay put.
        candidates.append((formula[:-1], 1))
    for candidate, protons in candidates:
        charge = ANION_CHARGES.get(candidate)
        if charge == -protons and _composition_of(candidate) is not None:
            return (candidate, charge, protons)
    return None


def _carbonate_parts(sp: Species) -> tuple[str, int, str, int] | None:
    """Split any carbonate or bicarbonate salt, whatever its cation."""
    parts = salt_parts(sp)
    if parts is None or parts[2] not in ("CO3", "HCO3"):
        return None
    return parts


def _hydroxide_parts(sp: Species) -> tuple[str, int, int] | None:
    counts = sp.counts
    if "O" not in counts or "H" not in counts:
        return None
    rest = {sym: n for sym, n in counts.items() if sym not in ("O", "H")}
    if len(rest) != 1:
        return None
    metal, n_metal = next(iter(rest.items()))
    if metal not in METALS:
        return None
    n_oh = counts["O"]
    if counts["H"] != n_oh or n_metal <= 0:
        return None
    states = [s for s in ELEMENTS[metal].states if s > 0]
    if not states:
        return None
    if n_oh % n_metal:
        return None
    # The hydroxide itself says how charged the cation must be: Fe(OH)3 -> +3.
    charge = int(n_oh // n_metal)
    if charge not in states:
        return None
    return (metal, charge, int(n_oh))


def _combine(cation: str, cation_charge: int, anion: str, anion_charge: int) -> str:
    return formula_from_ions(cation, cation_charge, anion, anion_charge)


def _key(text: str) -> str:
    """Match key for a formula as it is normally written (see ``match_key``)."""
    return match_key(parse_species(text))


_AMMONIA_KEY = _key("NH3")
_WATER_KEY = _key("H2O")
_HYDROGEN_KEY = _key("H2")


# Strong oxidisers that release a halogen from its acid.  Keyed by formula;
# each entry is (reduced metal, its state, E° couple, spectator cation).
_OXIDISERS: dict[str, tuple[str, int, str, str | None]] = {
    _key("KMnO4"): ("Mn", 2, "MnO4-/Mn2+", "K"),
    _key("NaMnO4"): ("Mn", 2, "MnO4-/Mn2+", "Na"),
    _key("K2Cr2O7"): ("Cr", 3, "Cr2O7^2-/Cr3+", "K"),
    _key("Na2Cr2O7"): ("Cr", 3, "Cr2O7^2-/Cr3+", "Na"),
    _key("MnO2"): ("Mn", 2, "MnO2/Mn2+", None),
    _key("PbO2"): ("Pb", 2, "PbO2/Pb2+", None),
}

_HALIDE_COUPLES = {"Cl": "Cl2/Cl-", "Br": "Br2/Br-", "I": "I2/I-"}
_HALOGEN_FORMULA = {"Cl": "Cl2", "Br": "Br2", "I": "I2"}

# Acidic (non-metal) oxides and the anion they form with a base.
_ACIDIC_OXIDES: dict[str, str] = {
    _key("CO2"): "CO3",
    _key("SO2"): "SO3",
    _key("SO3"): "SO4",
    _key("N2O5"): "NO3",
    _key("P4O10"): "PO4",
    _key("P2O5"): "PO4",
    _key("Cl2O7"): "ClO4",
}


def _oxide_parts(sp: Species) -> tuple[str, int] | None:
    """(metal, oxidation state) for a basic metal oxide such as Fe2O3."""
    counts = sp.counts
    if sp.charge != 0 or "O" not in counts or len(counts) != 2:
        return None
    if match_key(sp) in _OXIDISERS:
        return None  # MnO2/PbO2 are oxidisers, not plain basic oxides
    metal = next((sym for sym in counts if sym != "O"), None)
    if metal is None or metal not in METALS:
        return None
    oxygen = counts["O"]
    metal_atoms = counts[metal]
    if oxygen <= 0 or metal_atoms <= 0:
        return None
    charge = 2.0 * oxygen / metal_atoms
    if abs(charge - round(charge)) > 1e-9:
        return None
    state = int(round(charge))
    if state <= 0 or state > 4:
        return None  # Mn2O7, CrO3 and friends are acidic oxides, not basic ones
    if state not in [s for s in ELEMENTS[metal].states if s > 0]:
        return None
    return (metal, state)


def predict(reactants: list[Species], env=None) -> Prediction:
    """Predict products for a set of reactants."""
    input_keys = {match_key(sp) for sp in reactants}
    if not input_keys:
        raise PredictionError("no reactants given")

    best: _Entry | None = None
    for keys, entry in _ENTRY_TABLE:
        if keys == input_keys:
            best = entry
            break
    if best is not None:
        products = [parse_species(text) for text in best.products]
        return Prediction(products, "reaction database", [best.note])

    notes: list[str] = []
    products: list[Species] = []

    symbols = {sp: _element_symbol(sp) for sp in reactants}
    elements = [s for s in symbols.values() if s]

    # Combustion of any C/H/O compound in oxygen.
    if "O" in elements:
        fuel = next(
            (
                sp
                for sp in reactants
                if _element_symbol(sp) != "O"
                and "C" in sp.counts
                and set(sp.counts) <= {"C", "H", "O"}
            ),
            None,
        )
        if fuel is not None and set(fuel.counts) != {"C", "O"}:
            if env is not None and env.oxygen_fraction <= 0.001:
                raise PredictionError(
                    "this is a combustion reaction but the surroundings contain no oxygen"
                )
            products = [parse_species("CO2"), parse_species("H2O")]
            return Prediction(products, "combustion heuristic", ["assumes complete combustion"])

    # Metal + nonmetal -> binary ionic compound.
    metal_species = next((sp for sp in reactants if symbols[sp] in METALS and _is_element(sp)), None)
    nonmetal_species = next(
        (
            sp
            for sp in reactants
            if symbols[sp] in NONMETALS
            and _is_element(sp)
            and symbols[sp] not in ("H", "O")
        ),
        None,
    )
    if metal_species is not None and nonmetal_species is not None:
        metal = symbols[metal_species]
        nonmetal = symbols[nonmetal_species]
        cation_charge = ELEMENTS[metal].states[0] if ELEMENTS[metal].states else None
        anion_charge = _anion_charge(nonmetal)
        if cation_charge and anion_charge:
            formula = _combine(metal, cation_charge, nonmetal, anion_charge)
            return Prediction(
                [parse_species(formula)],
                "ionic compound heuristic",
                [f"{metal} is oxidised, {nonmetal} is reduced"],
            )

    # Metal + oxygen -> oxide (when oxygen present as a reactant).
    if "O" in elements and metal_species is not None:
        metal = symbols[metal_species]
        cation_charge = ELEMENTS[metal].states[0] if ELEMENTS[metal].states else None
        if cation_charge:
            formula = _combine(metal, cation_charge, "O", -2)
            return Prediction([parse_species(formula)], "oxide heuristic", notes)

    # Acid + base.
    acid = next((sp for sp in reactants if _acid_parts(sp)), None)
    base = next((sp for sp in reactants if _hydroxide_parts(sp)), None)
    carbonate = next((sp for sp in reactants if _carbonate_parts(sp) is not None), None)
    oxidiser = next((sp for sp in reactants if match_key(sp) in _OXIDISERS), None)
    oxide = next((sp for sp in reactants if _oxide_parts(sp) is not None), None)
    acidic_oxide = next((sp for sp in reactants if match_key(sp) in _ACIDIC_OXIDES), None)
    halide_acid: Species | None = None
    halide_anion = ""
    for sp in reactants:
        piece = _acid_parts(sp)
        if piece is not None and piece[0] in _HALIDE_COUPLES:
            halide_acid, halide_anion = sp, piece[0]
            break
    if acid and base:
        a_anion, a_charge, _ = _acid_parts(acid)
        b_cation, b_charge, _ = _hydroxide_parts(base)
        salt = _combine(b_cation, b_charge, a_anion, a_charge)
        return Prediction(
            [parse_species(salt), parse_species("H2O")],
            "acid-base neutralisation heuristic",
            ["salt plus water"],
        )
    if acid and any(match_key(sp) == _AMMONIA_KEY for sp in reactants):
        a_anion, a_charge, _ = _acid_parts(acid)
        salt = _combine("NH4", 1, a_anion, a_charge)
        return Prediction([parse_species(salt)], "ammonium salt heuristic", notes)

    # Acidified oxidiser: permanganate, dichromate, MnO2 or PbO2 in a halogen
    # acid release the halogen and drop to their reduced cation.
    if oxidiser is not None and halide_acid is not None and len(reactants) == 2:
        reduced_metal, reduced_state, couple, counter = _OXIDISERS[match_key(oxidiser)]
        products = []
        if counter is not None:
            products.append(
                parse_species(_combine(counter, 1, halide_anion, -1))
            )
        products.append(
            parse_species(_combine(reduced_metal, reduced_state, halide_anion, -1))
        )
        products.append(parse_species("H2O"))
        products.append(parse_species(_HALOGEN_FORMULA[halide_anion]))
        delta = cell_potential(couple, _HALIDE_COUPLES[halide_anion])
        if delta is None:
            note = f"{oxidiser.formula} oxidises {halide_anion}− to {_HALOGEN_FORMULA[halide_anion]}"
        elif delta >= 0.1:
            note = (
                f"{oxidiser.formula} oxidises {halide_anion}− readily "
                f"(E°cell = {delta:+.2f} V)"
            )
        elif delta >= 0.0:
            note = (
                f"marginally favourable (E°cell = {delta:+.2f} V); concentrated "
                "acid and warming drive it along"
            )
        else:
            note = (
                f"unfavourable at standard conditions (E°cell = {delta:+.2f} V); "
                f"only hot concentrated {halide_acid.formula} drives it"
            )
        return Prediction(
            products,
            "acidified oxidiser heuristic",
            [note, f"{reduced_metal} is reduced to its +{reduced_state} state"],
        )

    # Metal oxide + acid -> salt + water (basic oxides neutralise acids).
    if oxide is not None and acid is not None and len(reactants) == 2:
        metal, state = _oxide_parts(oxide)
        a_anion, a_charge, _ = _acid_parts(acid)
        salt = _combine(metal, state, a_anion, a_charge)
        return Prediction(
            [parse_species(salt), parse_species("H2O")],
            "metal oxide + acid heuristic",
            [f"{oxide.formula} is a basic oxide, so {acid.formula} gives the salt plus water"],
        )

    # Acidic oxide + base -> salt + water.
    if acidic_oxide is not None and base is not None and len(reactants) == 2:
        anion = _ACIDIC_OXIDES[match_key(acidic_oxide)]
        b_cation, b_charge, _ = _hydroxide_parts(base)
        salt = _combine(b_cation, b_charge, anion, ANION_CHARGES[anion])
        return Prediction(
            [parse_species(salt), parse_species("H2O")],
            "acidic oxide + base heuristic",
            [f"{acidic_oxide.formula} is an acidic oxide, so it forms {anion} salts"],
        )

    # Acid + sulfite / sulfide -> the gas plus the salt.
    if acid is not None and len(reactants) == 2:
        a_anion, a_charge, _ = _acid_parts(acid)
        other = next((sp for sp in reactants if sp is not acid), None)
        other_parts = salt_parts(other) if other is not None else None
        if other_parts is not None:
            cation, c_charge, anion, _ = other_parts
            salt = parse_species(_combine(cation, c_charge, a_anion, a_charge))
            if anion == "SO3":
                return Prediction(
                    [salt, parse_species("H2O"), parse_species("SO2")],
                    "acid + sulfite heuristic",
                    ["sulfur dioxide is released"],
                )
            if anion == "S":
                return Prediction(
                    [salt, parse_species("H2S")],
                    "acid + sulfide heuristic",
                    ["hydrogen sulfide is released (toxic, rotten-egg smell)"],
                )

    # Metal oxide + hydrogen -> metal + water, for oxides hydrogen can reduce.
    if oxide is not None and len(reactants) == 2 and "H" in elements:
        hydrogen = next((sp for sp in reactants if match_key(sp) == _HYDROGEN_KEY), None)
        metal, _state = _oxide_parts(oxide)
        if hydrogen is not None and metal in REACTIVITY_SERIES and "Fe" in REACTIVITY_SERIES:
            if REACTIVITY_SERIES.index(metal) >= REACTIVITY_SERIES.index("Fe"):
                return Prediction(
                    [parse_species(metal), parse_species("H2O")],
                    "oxide + hydrogen heuristic",
                    [
                        f"hydrogen reduces {oxide.formula} to {metal}; "
                        "more reactive metals keep their oxygen"
                    ],
                )

    # Acid + carbonate / bicarbonate: any carbonate, not a fixed list.
    if carbonate is not None:
        if acid is not None:
            parts = _carbonate_parts(carbonate)
            a_anion, a_charge, _ = _acid_parts(acid)
            if parts:
                cation, c_charge, _, _ = parts
                salt = _combine(cation, c_charge, a_anion, a_charge)
                return Prediction(
                    [parse_species(salt), parse_species("H2O"), parse_species("CO2")],
                    "acid + carbonate heuristic",
                    ["carbon dioxide is released"],
                )

    # Acid + metal (activity series).
    metal_element = next(
        (sp for sp in reactants if symbols[sp] in METALS and _is_element(sp)), None
    )
    if acid is not None and metal_element is not None:
        metal = symbols[metal_element]
        if metal in REACTIVITY_SERIES and "H" in REACTIVITY_SERIES:
            if REACTIVITY_SERIES.index(metal) > REACTIVITY_SERIES.index("H"):
                raise PredictionError(
                    f"{ELEMENTS[metal].name.lower()} sits below hydrogen in the activity "
                    "series, so it does not displace H+ from dilute acid"
                )
        a_anion, a_charge, protons = _acid_parts(acid)
        charge = ELEMENTS[metal].states[0] if ELEMENTS[metal].states else None
        if charge:
            salt = _combine(metal, charge, a_anion, a_charge)
            return Prediction(
                [parse_species(salt), parse_species("H2")],
                "metal + acid heuristic",
                ["hydrogen gas is released"],
            )

    # Metal oxide + water -> hydroxide.
    metal_oxide = next(
        (
            sp
            for sp in reactants
            if set(sp.counts) == {"O"} or (len(sp.counts) == 2 and "O" in sp.counts
                                           and any(e in METALS for e in sp.counts))
        ),
        None,
    )
    if metal_oxide is not None and any(match_key(sp) == _WATER_KEY for sp in reactants):
        metal = next((e for e in metal_oxide.counts if e in METALS), None)
        if metal and ELEMENTS[metal].states:
            charge = ELEMENTS[metal].states[0]
            formula = _combine(metal, charge, "OH", -1)
            return Prediction([parse_species(formula)], "oxide + water heuristic", notes)

    # Salt + salt / salt + hydroxide: double displacement.
    salts = [sp for sp in reactants if salt_parts(sp) is not None]
    if len(salts) >= 2:
        first, second = salts[0], salts[1]
        p1 = salt_parts(first)
        p2 = salt_parts(second)
        if p1 and p2:
            c1, c1q, a1, a1q = p1
            c2, c2q, a2, a2q = p2
            formula_a = _combine(c1, c1q, a2, a2q)
            formula_b = _combine(c2, c2q, a1, a1q)
            if not is_soluble(c1, a2) or not is_soluble(c2, a1):
                products = [parse_species(formula_a), parse_species(formula_b)]
                return Prediction(
                    products,
                    "double displacement heuristic",
                    ["one product is insoluble and precipitates"],
                )
            raise PredictionError(
                "both possible products are soluble, so no precipitate forms "
                "(all ions stay in solution)"
            )

    # Single displacement: metal + salt of a less reactive metal.
    if metal_element is not None:
        metal = symbols[metal_element]
        salt = next(
            (sp for sp in reactants if salt_parts(sp) is not None and sp is not metal_element),
            None,
        )
        if salt is not None:
            parts = salt_parts(salt)
            if parts:
                cation, c_charge, anion, a_charge = parts
                if cation in REACTIVITY_SERIES and metal in REACTIVITY_SERIES:
                    if REACTIVITY_SERIES.index(metal) < REACTIVITY_SERIES.index(cation):
                        charge = ELEMENTS[metal].states[0] if ELEMENTS[metal].states else None
                        if charge:
                            new_salt = _combine(metal, charge, anion, a_charge)
                            return Prediction(
                                [parse_species(new_salt), parse_species(cation)],
                                "single displacement heuristic",
                                [f"{metal} is more reactive than {cation}"],
                            )
                    else:
                        raise PredictionError(
                            f"{ELEMENTS[metal].name.lower()} is not reactive enough to "
                            f"displace {cation} from solution"
                        )

    # One reactant with no database entry: thermal decomposition.
    if len(reactants) == 1:
        single = reactants[0]
        parts = salt_parts(single)
        if parts is not None:
            cation, c_charge, anion, _ = parts
            if anion == "CO3" and cation in METALS and cation not in ALKALI:
                oxide_formula = _combine(cation, c_charge, "O", -2)
                return Prediction(
                    [parse_species(oxide_formula), parse_species("CO2")],
                    "carbonate decomposition heuristic",
                    ["needs strong heating; alkali carbonates are stable"],
                )
        hydroxide = _hydroxide_parts(single)
        if hydroxide is not None:
            metal, charge, _ = hydroxide
            if metal not in ALKALI:
                oxide_formula = _combine(metal, charge, "O", -2)
                return Prediction(
                    [parse_species(oxide_formula), parse_species("H2O")],
                    "hydroxide decomposition heuristic",
                    ["needs heating; the oxide and steam are left"],
                )

    raise PredictionError(
        "I do not know this reaction. Give the products explicitly with '->', "
        "for example: Fe + S -> FeS"
    )
