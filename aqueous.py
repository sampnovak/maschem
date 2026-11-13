"""Aqueous chemistry: pH from strong/weak acids, bases and salt hydrolysis."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .data import KA, KW_298, STRONG_ACIDS, STRONG_BASES, WEAK_BASES
from .elements import ALKALI, ALKALINE, ELEMENTS, METALS
from .formula import Species
from .surroundings import Surroundings

CONJUGATE_BASE: dict[str, tuple[str, int]] = {
    # anion formula -> (parent acid as listed in KA, which Ka index)
    "CH3COO": ("CH3COOH", 0),
    "C2H3O2": ("CH3COOH", 0),
    "HCOO": ("HCOOH", 0),
    "CO3": ("H2CO3", 1),
    "HCO3": ("H2CO3", 0),
    "CN": ("HCN", 0),
    "S": ("H2S", 1),
    "HS": ("H2S", 0),
    "F": ("HF", 0),
    "ClO": ("HClO", 0),
    "BrO": ("HBrO", 0),
    "NO2": ("HNO2", 0),
    "PO4": ("H3PO4", 2),
    "HPO4": ("H3PO4", 1),
    "SO3": ("H2SO3", 1),
}

_ACID_ANION_KA = {
    "HCO3": 4.7e-11,
    "HS": 1.3e-13,
    "HPO4": 4.2e-13,
    "H2PO4": 6.3e-8,
    "HSO3": 6.3e-8,
    "HSO4": 1.2e-2,
}


@dataclass
class PHResult:
    ph: float | None
    details: list[str] = field(default_factory=list)
    neutral_ph: float = 7.0
    temperature_k: float = 298.15
    reason: str | None = None


def _hydroxide_count(sp: Species) -> int:
    counts = sp.counts
    if "O" not in counts or "H" not in counts:
        return 0
    if all(symbol in ("O", "H") or symbol in METALS for symbol in counts):
        oxygen = counts.get("O", 0)
        if counts.get("H", 0) == oxygen and oxygen > 0:
            return int(oxygen)
    return 0


def _acid_protons(sp: Species) -> int:
    if sp.formula not in STRONG_ACIDS:
        return 0
    protons = int(sp.counts.get("H", 0))
    if sp.formula == "H2SO4":
        return 1  # second proton is weak, handled separately
    return protons


def _salt_hydrolysis(sp: Species) -> tuple[str, list[tuple[float, float]], list[tuple[float, float]]]:
    """Return (role, weak_acids, weak_bases) entries for a dissolved salt."""
    from .reactions import salt_parts

    parts = salt_parts(sp)
    if not parts:
        return ("spectator", [], [])
    cation, _, anion, _ = parts
    acids: list[tuple[float, float]] = []
    bases: list[tuple[float, float]] = []
    role = "neutral salt"

    ka_bh = KA.get(cation)
    if cation == "NH4":
        acids.append((5.6e-10, 1))
        role = "acidic salt (hydrolysis of NH4+)"
    elif ka_bh:
        acids.append((ka_bh[0], 1))
        role = "acidic salt"

    if anion in CONJUGATE_BASE:
        acid_formula, index = CONJUGATE_BASE[anion]
        ka_values = KA.get(acid_formula)
        if ka_values:
            ka = ka_values[min(index, len(ka_values) - 1)]
            bases.append((KW_298 / ka, 1))
            role = "basic salt (hydrolysis of a weak-acid anion)"
    elif anion in _ACID_ANION_KA:
        ka = _ACID_ANION_KA[anion]
        bases.append((KW_298 / ka, 1))
        role = "basic salt"

    return (role, acids, bases)


def analyse(
    species_moles: list[tuple[Species, float]],
    volume_l: float,
    env: Surroundings,
) -> PHResult:
    neutral = env.p_kw() / 2.0
    if not env.water:
        return PHResult(
            None,
            neutral_ph=neutral,
            temperature_k=env.temperature_k,
            reason="no liquid water in these surroundings, so pH is undefined",
        )
    if volume_l <= 0:
        return PHResult(None, reason="solution volume must be positive")

    kw = 10 ** (-env.p_kw())
    strong_acid = 0.0
    strong_base = 0.0
    weak_acids: list[tuple[float, float]] = []  # (concentration, Ka)
    weak_bases: list[tuple[float, float]] = []  # (concentration, Kb)
    details: list[str] = []

    for sp, moles in species_moles:
        if moles <= 0:
            continue
        concentration = moles / volume_l
        formula = sp.formula

        protons = _acid_protons(sp)
        if protons:
            strong_acid += concentration * protons
            details.append(f"{formula} is a strong acid ({protons} H+)")
            if formula == "H2SO4":
                weak_acids.append((concentration, 1.2e-2))
            continue

        if formula in STRONG_BASES:
            count = _hydroxide_count(sp) or 1
            strong_base += concentration * count
            details.append(f"{formula} is a strong base ({count} OH-)")
            continue

        if formula in WEAK_BASES:
            weak_bases.append((concentration, WEAK_BASES[formula]))
            details.append(f"{formula} is a weak base")
            continue

        ka_values = KA.get(formula)
        if ka_values and sp.formula not in STRONG_ACIDS:
            if sp.charge < 0:
                weak_bases.append((concentration, kw / ka_values[0]))
                details.append(f"{formula} is the conjugate base of a weak acid")
            else:
                weak_acids.append((concentration, ka_values[0]))
                details.append(f"{formula} is a weak acid (Ka = {ka_values[0]:.1e})")
            continue

        if sp.charge == 0 and len(sp.counts) >= 2:
            role, salt_acids, salt_bases = _salt_hydrolysis(sp)
            for ka in salt_acids:
                weak_acids.append((concentration, ka[0]))
            for kb in salt_bases:
                weak_bases.append((concentration, kb[0]))
            if salt_acids or salt_bases:
                details.append(f"{formula}: {role}")
            continue

    # Charge balance: [H+] + [B-H+] = [OH-] + [A-] + strong acid anions
    def f(h: float) -> float:
        oh = kw / h
        positive = h + strong_base
        for concentration, kb in weak_bases:
            ka_bh = kw / kb
            positive += concentration * h / (h + ka_bh)
        negative = oh + strong_acid
        for concentration, ka in weak_acids:
            negative += concentration * ka / (ka + h)
        return positive - negative

    lo, hi = -16.0, 2.0
    f_lo = f(10 ** lo)
    f_hi = f(10 ** hi)
    if f_lo * f_hi > 0:
        ph = neutral
    else:
        for _ in range(200):
            mid = (lo + hi) / 2
            value = f(10 ** mid)
            if abs(value) < 1e-14:
                break
            if value * f_lo <= 0:
                hi = mid
            else:
                lo = mid
                f_lo = value
        ph = -math.log10(10 ** ((lo + hi) / 2))

    if not details:
        details.append("no acid, base or hydrolysable salt present")

    return PHResult(
        ph=ph,
        details=details,
        neutral_ph=neutral,
        temperature_k=env.temperature_k,
    )
