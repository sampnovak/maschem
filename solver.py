"""Turn a typed reaction or mixture into a full analysis."""

from __future__ import annotations

import re
from dataclasses import dataclass, field, replace

from . import aqueous, redox, thermo
from .balance import BalanceError, balance, equation_text
from .data import AVOGADRO, KA, STRONG_ACIDS, STRONG_BASES, WEAK_BASES
from .elements import METALS
from .formula import FormulaError, Species, parse_species
from .reactions import PredictionError, is_soluble, predict, salt_parts
from .surroundings import Surroundings

_ARROW = re.compile(r"\s*(?:-->|->|=>|\u2192|\u27f6|\s=\s)\s*")

_UNIT_WORDS = [
    (r"\bkilograms?\b", "kg"),
    (r"\bmilligrams?\b", "mg"),
    (r"\bgrams?\b", "g"),
    (r"\bmillimoles?\b", "mmol"),
    (r"\bmoles?\b", "mol"),
]

_AMOUNT = re.compile(
    r"^\s*(?P<num>\d+(?:\.\d+)?|\.\d+)\s*"
    r"(?P<unit>kg|mg|g|mol|mmol)?\s*"
    r"(?:of\s+)?(?P<rest>.+?)\s*$",
    re.IGNORECASE,
)

_UNIT_TO_MOLES = {"mol": 1.0, "mmol": 1e-3}
_MASS_TO_GRAMS = {"kg": 1000.0, "g": 1.0, "mg": 1e-3}


@dataclass
class Component:
    species: Species
    moles: float | None = None
    given: str | None = None
    coeff: int = 1


@dataclass
class Analysis:
    input_text: str
    env: Surroundings
    volume_l: float
    error: str | None = None
    reactants: list[Component] = field(default_factory=list)
    products: list[Component] = field(default_factory=list)
    equation: str = ""
    source: str = ""
    source_notes: list[str] = field(default_factory=list)
    extent: float = 1.0
    per_mole: bool = True
    limiting: str | None = None
    dh_per: float | None = None
    dh_total: float | None = None
    missing_thermo: list[str] = field(default_factory=list)
    cp_total: float | None = None
    cp_estimated: list[str] = field(default_factory=list)
    t_initial: float = 298.15
    t_peak: float | None = None
    t_final: float = 298.15
    heat_capacity: float | None = None
    cooling_tau: float | None = None
    loss_coefficient: float | None = None
    temperature_curve: list[tuple[float, float]] = field(default_factory=list)
    ph_curve: list[tuple[float, float]] = field(default_factory=list)
    curve_window: float | None = None
    sample_step: float | None = None
    timeseries: list[dict] = field(default_factory=list)
    initial_ph: aqueous.PHResult | None = None
    final_ph: aqueous.PHResult | None = None
    redox: redox.RedoxResult | None = None
    product_rows: list[dict] = field(default_factory=list)
    unreacted: list[dict] = field(default_factory=list)
    element_totals: dict[str, dict] = field(default_factory=dict)
    colour_before: str = ""
    colour_after: str = ""
    mass_in: float = 0.0
    mass_out: float = 0.0
    peak_phases: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def _normalise_units(text: str) -> str:
    for pattern, replacement in _UNIT_WORDS:
        text = re.sub(pattern, replacement, text, flags=re.IGNORECASE)
    return text


def _parse_component(token: str) -> Component:
    token = _normalise_units(token.strip())
    match = _AMOUNT.match(token)
    if not match:
        return Component(species=parse_species(token))

    unit = match.group("unit")
    rest = match.group("rest").strip()
    species = parse_species(rest)
    if unit is None:
        # A bare number in front of a formula is a stoichiometric coefficient.
        return Component(species=species)

    number = float(match.group("num"))
    unit = unit.lower()

    if unit in _UNIT_TO_MOLES:
        moles = number * _UNIT_TO_MOLES[unit]
    else:
        grams = number * _MASS_TO_GRAMS.get(unit, 1.0)
        moles = grams / species.molar_mass
    return Component(species=species, moles=moles, given=f"{number:g} {unit} {species.formula}")


def _parse_side(side: str) -> list[Component]:
    side = re.sub(
        r"\s+(?:mixed with|together with|reacts with|plus|and)\s+",
        " + ",
        side,
        flags=re.IGNORECASE,
    )
    tokens = [t for t in re.split(r"\s+\+\s+", side.strip()) if t.strip()]
    return [_parse_component(t) for t in tokens]


def infer_state(sp: Species, env: Surroundings) -> str:
    """Guess the phase of a species in this environment."""
    if sp.state:
        return sp.state
    formula = sp.formula
    if formula == "H2O":
        if env.temperature_k < 273.15:
            return "s"
        if env.temperature_k >= env.water_boiling_point():
            return "g"
        return "l"

    if sp.charge and env.water:
        return "aq"  # a dissolved ion stays in solution

    if formula in STRONG_ACIDS or formula in STRONG_BASES:
        # Strong acids and bases are handled as solutions when water is around.
        if env.water:
            return "aq"
        return "l" if env.temperature_k > 273.15 else "s"

    if formula in KA or formula in WEAK_BASES:
        # Weak acids/bases are solutions too, unless they are volatile enough
        # to leave one: H2S and NH3 bubble out of water at room temperature.
        boiling = thermo.boiling_point_at(formula, env.pressure_pa)
        if boiling is not None and env.temperature_k >= boiling + 40.0 and not env.is_vacuum:
            return "g"
        if env.water:
            return "aq"
        return "l" if env.temperature_k > 273.15 else "s"

    parts = salt_parts(sp)
    if parts:
        cation, _, anion, _ = parts
        if env.water and is_soluble(cation, anion):
            return "aq"
        return "s"

    bp = thermo.boiling_point_at(formula, env.pressure_pa)
    if bp is not None and env.temperature_k >= bp and not env.is_vacuum:
        return "g"

    # Melting points are tabulated by formula (Br2, I2, ...), so check the
    # formula before falling back to the bare element symbol.
    mp = thermo.melting_point(formula)
    if mp is not None:
        return "s" if env.temperature_k < mp else "l"

    counts = sp.counts
    if len(counts) == 1:
        symbol = next(iter(counts))
        if symbol in METALS:
            return "s"
        element_mp = thermo.melting_point(symbol)
        if element_mp is not None and env.temperature_k < element_mp:
            return "s"

    return "s"


def _apply_states(items: list[Component], env: Surroundings) -> None:
    for component in items:
        component.species = component.species.with_state(infer_state(component.species, env))


def _colour_of(sp: Species, phase: str) -> str:
    from .data import COMPOUND_COLOURS, GAS_COLOURS, ION_COLOURS

    if phase == "g" and sp.formula in GAS_COLOURS:
        return GAS_COLOURS[sp.formula]
    if sp.charge and (sp.formula, sp.charge) in ION_COLOURS:
        return ION_COLOURS[(sp.formula, sp.charge)]
    if sp.formula in COMPOUND_COLOURS:
        return COMPOUND_COLOURS[sp.formula]
    if sp.formula in GAS_COLOURS and phase in ("g", "l"):
        return GAS_COLOURS[sp.formula]
    if phase == "aq":
        return "colourless"
    return "no colour data"


def _solution_colour(items: list[tuple[Species, float]]) -> str:
    from .data import ION_COLOURS

    found: list[str] = []
    for sp, moles in items:
        if moles <= 0 or sp.state not in ("aq", None):
            continue
        # A dissolved salt is coloured by its ions, so CuSO4 is a blue
        # solution even though the salt itself is white.
        candidates = [ION_COLOURS.get((sp.formula, sp.charge))]
        parts = salt_parts(sp)
        if parts:
            cation, c_charge, anion, a_charge = parts
            candidates.append(ION_COLOURS.get((cation, c_charge)))
            candidates.append(ION_COLOURS.get((anion, a_charge)))
        for colour in candidates:
            if colour and colour != "colourless" and colour not in found:
                found.append(colour)
    if not found:
        return "colourless solution"
    return "solution: " + ", ".join(found)


def _limiting(reactants: list[Component]) -> tuple[float, str | None]:
    available = [
        (c.moles / c.coeff, c.species.display)
        for c in reactants
        if c.moles is not None and c.coeff > 0
    ]
    if not available:
        return 1.0, None
    extent, name = min(available, key=lambda pair: pair[0])
    return extent, name


def _element_totals(
    products: list[Component], unreacted: list[dict]
) -> dict[str, dict]:
    from .elements import ELEMENTS

    totals: dict[str, dict] = {}

    def add(species: Species, moles: float) -> None:
        for symbol, count in species.counts.items():
            entry = totals.setdefault(symbol, {"moles": 0.0, "atoms": 0.0, "mass": 0.0})
            amount = count * moles
            entry["moles"] += amount
            entry["atoms"] += amount * AVOGADRO
            entry["mass"] += amount * ELEMENTS[symbol].mass

    for component in products:
        add(component.species, component.moles or 0.0)
    for row in unreacted:
        add(row["species"], row["moles"])
    return totals


def time_series(
    analysis: Analysis, step: float | None = None, end: float | None = None
) -> list[dict]:
    """Temperature and pH after the reaction, at any time resolution.

    The reaction itself is treated as instantaneous: the mixture starts at its
    adiabatic peak (or dip) and relaxes as T(t) = T_env + (T_peak - T_env)
    exp(-t/tau).  pH is recomputed at every one of those temperatures, so it
    moves because Kw and the acid equilibria are temperature dependent.
    """
    if analysis.t_peak is None or not analysis.cooling_tau:
        return []
    tau = analysis.cooling_tau
    ambient = analysis.env.temperature_k
    window = analysis.curve_window or thermo.curve_window(tau)
    rate = float(step if step else (analysis.sample_step or thermo.curve_step(tau)))
    rate = max(rate, 1e-6)
    finish = float(window if end is None else max(end, 0.0))
    if finish / rate > 400:
        rate = finish / 400.0
    mixture = [
        (row["species"], row["moles"])
        for row in analysis.product_rows
        if row["species"].state == "aq"
    ] + [
        (row["species"], row["moles"])
        for row in analysis.unreacted
        if row["species"].state == "aq"
    ]
    boiling = analysis.env.water_boiling_point()
    rows: list[dict] = []
    index = 0
    while index * rate <= finish + 1e-9:
        moment = index * rate
        temperature = thermo.cooling_temperature(analysis.t_peak, ambient, tau, moment)
        ph = None
        if mixture and analysis.env.water and temperature < boiling:
            try:
                ph = aqueous.analyse(
                    mixture,
                    analysis.volume_l,
                    replace(analysis.env, temperature_k=temperature),
                ).ph
            except Exception:  # a pH model that cannot handle one point
                ph = None
        rows.append(
            {
                "t": round(moment, 6),
                "T": temperature,
                "C": temperature - 273.15,
                "ph": ph,
            }
        )
        index += 1
    return rows


def analyse(
    text: str,
    env: Surroundings,
    volume_l: float = 1.0,
    loss: float = thermo.LOSS_COEFFICIENT,
) -> Analysis:
    """Analyse one typed reaction.

    ``loss`` is the heat-loss coefficient in W/K used for the cooling model
    (0.05 for a dewar, 0.5 for an open beaker, 5 for something fan-cooled).
    """
    result = Analysis(input_text=text, env=env, volume_l=volume_l)
    try:
        _analyse(text, env, volume_l, result, loss)
    except (FormulaError, BalanceError, PredictionError) as exc:
        result.error = str(exc)
    except ZeroDivisionError:
        result.error = "cannot compute: a reactant has zero molar mass"
    except Exception as exc:  # never let one odd input kill the shell
        result.error = f"could not analyse this input: {exc}"
    return result


def _analyse(
    text: str, env: Surroundings, volume_l: float, result: Analysis, loss: float
) -> None:
    text = text.strip()
    if not text:
        raise FormulaError("nothing to solve")

    parts = _ARROW.split(text, maxsplit=1)
    left_text = parts[0]
    right_text = parts[1] if len(parts) > 1 else ""

    reactants = _parse_side(left_text)
    if not reactants:
        raise FormulaError("no reactants recognised")

    if right_text.strip():
        products = _parse_side(right_text)
        source = "products given"
        source_notes: list[str] = []
        result.source = source
    else:
        prediction = predict([c.species for c in reactants], env)
        products = [Component(species=sp) for sp in prediction.products]
        source_notes = prediction.notes
        result.source = prediction.source
    result.source_notes = source_notes

    rcoeffs, pcoeffs, balance_notes = balance(
        [c.species for c in reactants], [c.species for c in products]
    )

    keep_r = [i for i, c in enumerate(rcoeffs) if c != 0]
    keep_p = [i for i, c in enumerate(pcoeffs) if c != 0]
    dropped = [
        reactants[i].species.display
        for i, c in enumerate(rcoeffs)
        if c == 0
    ] + [products[i].species.display for i, c in enumerate(pcoeffs) if c == 0]
    reactants = [reactants[i] for i in keep_r]
    rcoeffs = [rcoeffs[i] for i in keep_r]
    products = [products[i] for i in keep_p]
    pcoeffs = [pcoeffs[i] for i in keep_p]

    for component, coeff in zip(reactants, rcoeffs):
        component.coeff = coeff
    for component, coeff in zip(products, pcoeffs):
        component.coeff = coeff

    _apply_states(reactants, env)
    _apply_states(products, env)

    extent, limiting = _limiting(reactants)
    result.extent = extent
    result.limiting = limiting
    result.per_mole = all(c.moles is None for c in reactants)

    # Fill in moles for species without supplied amounts.
    for component in reactants:
        if component.moles is None:
            component.moles = component.coeff * extent
    for component in products:
        component.moles = component.coeff * extent

    result.reactants = reactants
    result.products = products
    result.equation = equation_text(
        [c.species for c in reactants], rcoeffs, [c.species for c in products], pcoeffs
    )

    if balance_notes:
        result.notes.extend(balance_notes)
    if dropped:
        result.notes.append("not part of the reaction: " + ", ".join(dropped))
    if result.per_mole:
        result.notes.append(
            "no amounts were given, so results are per mole of the balanced reaction"
        )
    if limiting:
        result.notes.append(f"limiting reagent: {limiting}")

    effective_volume = volume_l
    water_reactant = next(
        (c for c in reactants if c.species.formula == "H2O" and c.moles), None
    )
    if water_reactant is not None and water_reactant.moles:
        effective_volume = max(effective_volume, water_reactant.moles * 0.018)
    result.volume_l = effective_volume

    # Thermochemistry.
    dh, missing, estimated = thermo.enthalpy_of_reaction(
        [c.species for c in reactants], rcoeffs, [c.species for c in products], pcoeffs
    )
    result.dh_per = dh
    result.missing_thermo = missing
    aqueous_reaction = any(c.species.state == "aq" for c in reactants + products)
    if dh is not None:
        result.dh_total = dh * extent
        if aqueous_reaction and env.water:
            # The solution, not the stoichiometric product mixture, absorbs the heat.
            result.t_peak = env.temperature_k + (
                -result.dh_total * 1000.0
            ) / (4180.0 * effective_volume)
            result.notes.append(
                f"temperature is for {effective_volume:g} L of solution absorbing the heat "
                "(water heat capacity 4.18 J/g K)"
            )
        else:
            cp_total, cp_estimated = thermo.products_heat_capacity(
                [c.species for c in products], pcoeffs
            )
            result.cp_total = cp_total
            result.cp_estimated = cp_estimated
            result.t_peak = thermo.peak_temperature(
                dh, [c.species for c in products], pcoeffs, env.temperature_k
            )

    result.t_initial = env.temperature_k
    result.t_final = env.temperature_k

    # Temperature as it evolves after the reaction has dumped its heat in:
    # the mixture starts at its adiabatic peak and cools towards the
    # surroundings with the heat capacity it actually has.
    if result.t_peak is not None:
        if aqueous_reaction and env.water:
            heat_capacity = 4180.0 * effective_volume
        else:
            cp_total = result.cp_total
            if cp_total is None:
                cp_total, _ = thermo.products_heat_capacity(
                    [c.species for c in products], pcoeffs
                )
            heat_capacity = cp_total * extent
        result.heat_capacity = max(heat_capacity, 1e-9)
        result.loss_coefficient = loss
        result.cooling_tau = thermo.cooling_time_constant(result.heat_capacity, loss)
        # The automatic window follows the time constant, so a slow vessel
        # still shows its dip and recovery instead of 20 identical rows.
        result.curve_window = thermo.curve_window(result.cooling_tau)
        result.sample_step = thermo.curve_step(result.cooling_tau)

    # Redox.
    result.redox = redox.analyse(
        [c.species for c in reactants],
        rcoeffs,
        [c.species for c in products],
        pcoeffs,
    )

    # Product table.
    for component in products:
        sp = component.species
        moles = component.moles or 0.0
        result.product_rows.append(
            {
                "species": sp,
                "moles": moles,
                "mass": moles * sp.molar_mass,
                "atoms": moles * sp.atom_count * AVOGADRO,
            }
        )

    for component in reactants:
        if component.moles is None:
            continue
        if result.per_mole:
            continue
        left = component.moles - component.coeff * extent
        if left > 1e-12:
            result.unreacted.append(
                {
                    "species": component.species,
                    "moles": left,
                    "mass": left * component.species.molar_mass,
                }
            )

    result.element_totals = _element_totals(products, result.unreacted)

    result.mass_in = sum((c.moles or 0.0) * c.species.molar_mass for c in reactants)
    result.mass_out = sum(row["mass"] for row in result.product_rows) + sum(
        row["mass"] for row in result.unreacted
    )

    # Aqueous chemistry.
    if env.water:
        before = [
            (c.species, c.moles or 0.0) for c in reactants if c.species.state == "aq"
        ]
        after = [
            (c.species, c.moles or 0.0) for c in products if c.species.state == "aq"
        ]
        after += [
            (row["species"], row["moles"])
            for row in result.unreacted
            if row["species"].state == "aq"
        ]
        if before:
            result.initial_ph = aqueous.analyse(before, effective_volume, env)
        if after:
            result.final_ph = aqueous.analyse(after, effective_volume, env)
        result.colour_before = _solution_colour(before)
        result.colour_after = _solution_colour(after)
    else:
        result.colour_before = "not applicable (no solvent)"
        result.colour_after = "not applicable (no solvent)"

    if not env.water:
        result.notes.append("no liquid water: pH is not defined and ionic products stay solid")

    # Time series: temperature and pH on the same timetable.
    if result.t_peak is not None and result.cooling_tau:
        result.timeseries = time_series(result)
        result.temperature_curve = [(row["t"], row["T"]) for row in result.timeseries]
        result.ph_curve = [
            (row["t"], row["ph"]) for row in result.timeseries if row["ph"] is not None
        ]

    # Phases at the peak temperature.
    if result.t_peak is not None and not aqueous_reaction:
        for component in products:
            phase = thermo.phase_at(component.species, result.t_peak, env.pressure_pa)
            result.peak_phases.append(
                f"{component.species.display} \u2192 {phase} at {result.t_peak:.0f} K"
            )

    # Warnings.
    if any(c.species.formula == "O2" for c in reactants) and env.oxygen_fraction < 1e-3:
        result.warnings.append(
            "this reaction consumes oxygen but these surroundings contain none "
            "(you would need to supply an oxidiser)"
        )
    if env.is_vacuum:
        gases = [c.species.formula for c in products if c.species.state == "g"]
        if gases:
            result.warnings.append(
                "vacuum: gaseous products (" + ", ".join(gases) + ") escape immediately"
            )
        result.warnings.append(
            "vacuum: no convective cooling, heat is lost only by radiation"
        )
    if dh is not None and dh > 0 and not aqueous_reaction:
        result.warnings.append(
            "endothermic: this reaction needs continuous heating; it will not "
            "keep itself going"
        )
    if result.t_peak is None and dh is not None and dh < 0:
        result.warnings.append(
            "the reaction releases too little heat to keep itself going at these "
            "surroundings; the adiabatic temperature falls below absolute zero"
        )
    if result.t_peak is not None and result.t_peak > 3000:
        result.warnings.append(
            "the adiabatic temperature is an upper bound: dissociation and heat loss "
            "make the real temperature lower"
        )
    if result.dh_per is not None and result.dh_per < -2000:
        result.warnings.append("very exothermic; treat as explosive/hazardous")
    if result.missing_thermo:
        result.warnings.append(
            "no formation enthalpy for: " + ", ".join(sorted(set(result.missing_thermo)))
        )

    if result.redox is not None and not result.redox.happened:
        result.notes.append("no element changes oxidation state: not a redox reaction")
