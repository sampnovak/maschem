"""Chemical formula parsing for MasChem.

Handles element counts, nested brackets, hydrates (``CuSO4·5H2O``),
ionic charges (``SO4^2-``, ``Fe3+``) and phase labels (``(aq)``).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace

from .elements import ELEMENTS

STATE_RE = re.compile(r"\((aq|s|l|g)\)\s*$", re.IGNORECASE)
CARET_CHARGE_RE = re.compile(r"\^\s*(\d*)\s*([+-])\)?$")
BARE_CHARGE_RE = re.compile(r"(\d*)([+-])$")
LEADING_NUMBER_RE = re.compile(r"^\d+(?=[A-Z(\[])")
SINGLE_ELEMENT_RE = re.compile(r"^[A-Z][a-z]?$")


class FormulaError(ValueError):
    pass


# Common names people type instead of formulas.
NAME_ALIASES = {
    "water": "H2O",
    "ice": "H2O",
    "steam": "H2O",
    "hydrogen peroxide": "H2O2",
    "table salt": "NaCl",
    "salt": "NaCl",
    "iron": "Fe",
    "copper": "Cu",
    "zinc": "Zn",
    "aluminium": "Al",
    "aluminum": "Al",
    "magnesium": "Mg",
    "sodium": "Na",
    "potassium": "K",
    "calcium": "Ca",
    "carbon": "C",
    "sulfur": "S",
    "sulphur": "S",
    "oxygen": "O2",
    "hydrogen": "H2",
    "nitrogen": "N2",
    "chlorine": "Cl2",
    "methane": "CH4",
    "ethane": "C2H6",
    "propane": "C3H8",
    "butane": "C4H10",
    "ethanol": "C2H5OH",
    "methanol": "CH3OH",
    "benzene": "C6H6",
    "glucose": "C6H12O6",
    "acetic acid": "CH3COOH",
    "vinegar": "CH3COOH",
    "ammonia": "NH3",
    "carbon dioxide": "CO2",
    "carbon monoxide": "CO",
    "limestone": "CaCO3",
    "quicklime": "CaO",
    "slaked lime": "Ca(OH)2",
    "baking soda": "NaHCO3",
    "washing soda": "Na2CO3",
    "caustic soda": "NaOH",
    "lye": "NaOH",
    "rust": "Fe2O3",
    "quartz": "SiO2",
    "sand": "SiO2",
    "bleach": "NaClO",
    "nitric acid": "HNO3",
    "sulfuric acid": "H2SO4",
    "sulphuric acid": "H2SO4",
    "hydrochloric acid": "HCl",
    "phosphoric acid": "H3PO4",
    "potassium permanganate": "KMnO4",
    "copper sulfate": "CuSO4",
    "copper sulphate": "CuSO4",
}

_HYDRATE_SPLIT = re.compile(r"[\u00b7*]|\.(?=\d)")


@dataclass(frozen=True)
class Species:
    """A chemical species: formula + charge + optional phase."""

    formula: str
    composition: tuple[tuple[str, float], ...]
    charge: int = 0
    state: str | None = None

    @property
    def counts(self) -> dict[str, float]:
        return dict(self.composition)

    @property
    def molar_mass(self) -> float:
        total = 0.0
        for sym, n in self.composition:
            if sym not in ELEMENTS:
                raise FormulaError(f"unknown element {sym!r} in {self.formula!r}")
            total += ELEMENTS[sym].mass * n
        return total

    @property
    def atom_count(self) -> float:
        return sum(n for _, n in self.composition)

    @property
    def elements(self) -> list[str]:
        return [s for s, _ in self.composition]

    @property
    def charge_text(self) -> str:
        if self.charge == 0:
            return ""
        sign = "+" if self.charge > 0 else "\u2212"
        mag = abs(self.charge)
        return sign if mag == 1 else f"{mag}{sign}"

    @property
    def display(self) -> str:
        text = self.formula + self.charge_text
        if self.state:
            text += f"({self.state})"
        return text

    @property
    def state_text(self) -> str:
        return self.state or "?"

    def with_state(self, state: str | None) -> "Species":
        return replace(self, state=state)

    def __str__(self) -> str:  # pragma: no cover - convenience
        return self.display


def _parse_segment(segment: str) -> dict[str, float]:
    counts: dict[str, float] = {}
    pos = 0
    for match in re.finditer(r"([A-Z][a-z]?)(\d*\.?\d*)", segment):
        if match.start() != pos:
            raise FormulaError(f"cannot parse {segment!r} near {segment[pos:]!r}")
        symbol = match.group(1)
        if symbol not in ELEMENTS:
            raise FormulaError(f"unknown element {symbol!r}")
        number = match.group(2)
        amount = float(number) if number else 1.0
        counts[symbol] = counts.get(symbol, 0.0) + amount
        pos = match.end()
    if pos != len(segment):
        raise FormulaError(f"cannot parse {segment!r} near {segment[pos:]!r}")
    if not counts:
        raise FormulaError(f"empty formula {segment!r}")
    return counts


def _expand_groups(text: str) -> str:
    while "(" in text or "[" in text:
        match = re.search(r"[\(\[]([^\(\)\[\]]*)[\)\]](\d*\.?\d*)", text)
        if not match:
            raise FormulaError(f"unbalanced brackets in {text!r}")
        inner = match.group(1)
        multiplier = float(match.group(2) or 1)
        counts = _parse_segment(inner)
        pieces = []
        for symbol, amount in counts.items():
            total = amount * multiplier
            if abs(total - round(total)) < 1e-9:
                total = int(round(total))
            pieces.append(symbol if total == 1 else f"{symbol}{total}")
        text = text[: match.start()] + "".join(pieces) + text[match.end() :]
    return text


def parse_formula(text: str) -> dict[str, float]:
    """Parse a formula (no charge/state) into an element -> count mapping."""
    text = text.strip().replace(" ", "")
    if not text:
        raise FormulaError("empty formula")
    parts = _HYDRATE_SPLIT.split(text)
    combined: dict[str, float] = {}
    for part in parts:
        if not part:
            continue
        match = re.match(r"^(\d+(?:\.\d+)?)(.*)$", part)
        if match and match.group(2) and ("(" in match.group(2) or match.group(2)[:1].isupper()):
            multiplier = float(match.group(1))
            body = match.group(2)
        else:
            multiplier = 1.0
            body = part
        counts = _parse_segment(_expand_groups(body))
        for symbol, amount in counts.items():
            combined[symbol] = combined.get(symbol, 0.0) + amount * multiplier
    return combined


def parse_species(text: str, state: str | None = None) -> Species:
    raw = text.strip()
    if not raw:
        raise FormulaError("empty species")

    alias = NAME_ALIASES.get(raw.lower())
    if alias:
        raw = alias

    match = STATE_RE.search(raw)
    if match:
        state = state or match.group(1).lower()
        raw = raw[: match.start()].strip()

    charge = 0
    match = CARET_CHARGE_RE.search(raw)
    if match:
        # An explicit caret is unambiguous: SO4^2-, Fe(CN)6^3-, NH4^+.
        magnitude = int(match.group(1) or 1)
        charge = magnitude if match.group(2) == "+" else -magnitude
        raw = raw[: match.start()].strip()
    else:
        match = BARE_CHARGE_RE.search(raw)
        if match:
            digits = match.group(1)
            body = raw[: match.start()].strip()
            sign = 1 if match.group(2) == "+" else -1
            if digits and SINGLE_ELEMENT_RE.match(body):
                # Fe3+, Cu2+, O2-: the digits are the charge on one atom.
                charge = sign * int(digits)
                raw = body
            else:
                # NH4+, H3O+, NO3-: any digits are subscripts, charge is 1.
                charge = sign
                raw = raw[:-1].strip()

    raw = LEADING_NUMBER_RE.sub("", raw, count=1)
    if not raw:
        raise FormulaError(f"cannot read species {text!r}")

    composition = parse_formula(raw)
    frozen = tuple(sorted(composition.items()))
    return Species(formula=raw, composition=frozen, charge=charge, state=state)


def hill_formula(counts: dict[str, float]) -> str:
    """Canonical formula for matching (Hill order)."""
    if not counts:
        return ""

    def token(symbol: str, amount: float) -> str:
        if abs(amount - 1) < 1e-9:
            return symbol
        if abs(amount - round(amount)) < 1e-9:
            return f"{symbol}{int(round(amount))}"
        return f"{symbol}{amount:g}"

    if "C" in counts:
        order = ["C"]
        if "H" in counts:
            order.append("H")
        order += sorted(s for s in counts if s not in ("C", "H"))
    else:
        order = sorted(counts)
    return "".join(token(s, counts[s]) for s in order)


def match_key(sp: Species) -> str:
    """Key used to look reactions up: element symbol or Hill formula."""
    counts = sp.counts
    if len(counts) == 1:
        return next(iter(counts))
    return hill_formula(counts)


def formula_from_ions(cation: str, cation_charge: int, anion: str, anion_charge: int) -> str:
    """Build a neutral ionic formula from ion formulas and charges."""
    from math import gcd

    a = abs(anion_charge)
    b = abs(cation_charge)
    g = gcd(a, b) or 1
    n_cat = a // g
    n_an = b // g
    return _scale_formula(cation, n_cat) + _scale_formula(anion, n_an)


def _scale_formula(formula: str, mult: int) -> str:
    if mult == 1:
        return formula
    counts = parse_formula(formula)
    if len(counts) == 1:
        symbol, amount = next(iter(counts.items()))
        total = amount * mult
        if abs(total - round(total)) < 1e-9:
            total = int(round(total))
        return symbol if total == 1 else f"{symbol}{total}"
    # A polyatomic ion keeps its identity: Al2(SO4)3, Ca(OH)2, (NH4)2SO4.
    return f"({formula}){mult}"
