"""Thermochemistry: reaction enthalpy, adiabatic temperature, phase changes."""

from __future__ import annotations

import math

from .data import BP, CP, DHF, MP, DH_VAP_WATER, R
from .formula import Species

STANDARD_STATE = {"H2", "N2", "O2", "F2", "Cl2", "Br2", "I2", "C", "S", "P4"}

# Heat lost to the room, in watts per kelvin, for a small open vessel cooling
# by natural convection.  Change it to model a flask, a vacuum dewar, etc.
LOSS_COEFFICIENT = 0.5


def cooling_time_constant(
    total_heat_capacity: float, loss: float = LOSS_COEFFICIENT
) -> float:
    """Newtonian cooling time constant, seconds: C / (loss coefficient)."""
    if total_heat_capacity <= 0 or loss <= 0:
        return 0.0
    return total_heat_capacity / loss


def cooling_temperature(
    peak_k: float, ambient_k: float, tau_s: float, time_s: float
) -> float:
    """Temperature ``time_s`` after the reaction dumped its heat in.

    The reaction is treated as fast compared with cooling, so the mixture
    starts at its adiabatic peak and then relaxes towards the surroundings:
    T(t) = T_env + (T_peak - T_env) * exp(-t / tau).
    """
    if tau_s <= 0:
        return ambient_k
    return ambient_k + (peak_k - ambient_k) * math.exp(-time_s / tau_s)


def cooling_curve(
    peak_k: float,
    ambient_k: float,
    tau_s: float,
    times: list[float],
) -> list[tuple[float, float]]:
    """(time, temperature) pairs for a list of sample times."""
    return [(t, cooling_temperature(peak_k, ambient_k, tau_s, t)) for t in times]


# Sampling steps that read nicely on a terminal, in seconds.
NICE_STEPS = (
    0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 30.0, 60.0, 120.0,
    300.0, 600.0, 900.0, 1800.0, 3600.0, 7200.0, 14400.0, 28800.0,
)
CURVE_FACTOR = 4.0  # four time constants is 98% of the way back
CURVE_ROWS = 24


def curve_window(tau_s: float | None, factor: float = CURVE_FACTOR) -> float:
    """How long a cooling curve is worth watching: four time constants.

    A one-second window would miss everything in a litre of water, which
    takes hours to cool; the window therefore follows the time constant.
    """
    if not tau_s or tau_s <= 0:
        return 0.0
    return tau_s * factor


def curve_step(
    tau_s: float | None, max_rows: int = CURVE_ROWS, factor: float = CURVE_FACTOR
) -> float:
    """Coarsest nice step that still shows the whole curve in ``max_rows``."""
    window = curve_window(tau_s, factor)
    if window <= 0:
        return 1.0
    for step in NICE_STEPS:
        if window / step <= max_rows:
            return step
    return NICE_STEPS[-1]


def curve_times(
    tau_s: float | None, factor: float = CURVE_FACTOR, max_rows: int = CURVE_ROWS
) -> list[float]:
    """Sample times from 0 to four tau, end point included."""
    window = curve_window(tau_s, factor)
    if window <= 0:
        return [0.0]
    step = curve_step(tau_s, max_rows, factor)
    times: list[float] = []
    index = 0
    while index * step <= window + 1e-9:
        times.append(round(index * step, 6))
        index += 1
    return times


def delta_hf(sp: Species) -> float | None:
    keys = []
    if sp.state:
        keys.append(f"{sp.formula}({sp.state})")
    keys.append(sp.formula)
    for key in keys:
        if key in DHF:
            return DHF[key]
    counts = sp.counts
    if len(counts) == 1:
        symbol = next(iter(counts))
        if sp.formula in STANDARD_STATE or sp.formula == symbol:
            return 0.0
    return None


def heat_capacity(sp: Species) -> tuple[float, bool]:
    """Cp in J/(mol K); second value says whether it was estimated."""
    keys = []
    if sp.state:
        keys.append(f"{sp.formula}({sp.state})")
    keys.append(sp.formula)
    for key in keys:
        if key in CP:
            return CP[key], False
    atoms = max(sp.atom_count, 1.0)
    state = sp.state or ("g" if sp.formula in BP and sp.formula in ("H2", "O2", "N2") else "")
    if state == "g":
        return 25.0 + 4.0 * atoms, True
    if state == "l":
        return 22.0 * atoms, True
    return 25.0 * atoms, True


def enthalpy_of_reaction(
    reactants: list[Species],
    rcoeffs: list[int],
    products: list[Species],
    pcoeffs: list[int],
) -> tuple[float | None, list[str], list[str]]:
    missing: list[str] = []
    total = 0.0
    for sp, coeff in zip(reactants, rcoeffs):
        if coeff == 0:
            continue
        value = delta_hf(sp)
        if value is None:
            missing.append(sp.display)
        else:
            total -= value * coeff
    for sp, coeff in zip(products, pcoeffs):
        if coeff == 0:
            continue
        value = delta_hf(sp)
        if value is None:
            missing.append(sp.display)
        else:
            total += value * coeff
    if missing:
        return None, missing, []
    estimated = []
    for sp in products:
        _, guess = heat_capacity(sp)
        if guess:
            estimated.append(sp.display)
    return total, [], estimated


def products_heat_capacity(
    products: list[Species], pcoeffs: list[int]
) -> tuple[float, list[str]]:
    total = 0.0
    estimated: list[str] = []
    for sp, coeff in zip(products, pcoeffs):
        if coeff == 0:
            continue
        cp, guess = heat_capacity(sp)
        total += cp * coeff
        if guess:
            estimated.append(sp.display)
    return total, estimated


def peak_temperature(
    dh_per_kj: float | None,
    products: list[Species],
    pcoeffs: list[int],
    start_k: float,
) -> float | None:
    """Adiabatic temperature for one mole of reaction.

    Uses a simple high-temperature heat-capacity correction and switches a
    product to its gas heat capacity above its boiling point.  Dissociation
    and heat loss are ignored, so this is an upper bound.
    """
    if dh_per_kj is None:
        return None

    def cp_at(temperature: float) -> float:
        total = 0.0
        for sp, coeff in zip(products, pcoeffs):
            if coeff == 0:
                continue
            cp, _ = heat_capacity(sp)
            bp = boiling_point_at(sp.formula, 101325.0)
            if bp is not None and temperature > bp:
                total += coeff * cp * max(1.0, (temperature / 298.15) ** 0.25)
            else:
                total += coeff * cp
        return total

    cp = cp_at(start_k)
    if cp <= 0:
        return None
    estimate = start_k - (dh_per_kj * 1000.0) / cp
    if estimate <= 1.0:
        return None
    for _ in range(60):
        cp = cp_at(estimate)
        if cp <= 0:
            return None
        updated = start_k - (dh_per_kj * 1000.0) / cp
        if updated <= 1.0:
            return None
        if abs(updated - estimate) < 0.5:
            estimate = updated
            break
        estimate = 0.5 * estimate + 0.5 * updated
    return estimate


def boiling_point_at(formula: str, pressure_pa: float) -> float | None:
    base = BP.get(formula)
    if base is None:
        return None
    if pressure_pa <= 0:
        return None
    if abs(pressure_pa - 101325.0) < 1.0:
        return base
    dh = DH_VAP_WATER if formula == "H2O" else 88.0 * base  # Trouton's rule
    ln = math.log(pressure_pa / 101325.0)
    inv = 1.0 / base - (R * ln) / dh
    if inv <= 0:
        return 0.0
    return 1.0 / inv


def melting_point(formula: str) -> float | None:
    return MP.get(formula)


def phase_at(sp: Species, temperature_k: float, pressure_pa: float) -> str:
    """Best guess of the phase at a temperature/pressure."""
    formula = sp.formula
    mp = melting_point(formula)
    bp = boiling_point_at(formula, pressure_pa)
    if bp is not None and temperature_k >= bp:
        return "g"
    if pressure_pa <= 0:
        if mp is not None:
            return "s" if temperature_k < mp else "g"
        return sp.state or "s"
    if mp is not None:
        return "s" if temperature_k < mp else "l"
    return sp.state or "s"
