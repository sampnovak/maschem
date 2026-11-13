"""Surrounding environments for MasChem (E1 Earth, S1 space, customs)."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from .data import BP, DH_ION_WATER, DH_VAP_WATER, KW_298, R

PRESETS: dict[str, "Surroundings"] = {}


@dataclass(frozen=True)
class Surroundings:
    code: str
    name: str
    temperature_k: float
    pressure_pa: float
    atmosphere: tuple[tuple[str, float], ...]  # mole fractions
    water: bool
    notes: tuple[str, ...] = ()

    # -- derived ---------------------------------------------------------
    @property
    def temperature_c(self) -> float:
        return self.temperature_k - 273.15

    @property
    def pressure_atm(self) -> float:
        return self.pressure_pa / 101325.0

    @property
    def is_vacuum(self) -> bool:
        return self.pressure_pa < 1.0

    @property
    def oxygen_fraction(self) -> float:
        return dict(self.atmosphere).get("O2", 0.0)

    def p_kw(self) -> float:
        """pKw at this temperature (van't Hoff, water ionisation)."""
        t = max(self.temperature_k, 1.0)
        ln_ratio = -(DH_ION_WATER / R) * (1.0 / t - 1.0 / 298.15)
        kw = KW_298 * pow(2.718281828459045, ln_ratio)
        return -__import__("math").log10(kw)

    def water_boiling_point(self) -> float:
        """Boiling point of water at this pressure (Clausius-Clapeyron)."""
        if self.pressure_pa <= 0:
            return 0.0
        if abs(self.pressure_pa - 101325.0) < 1e-6:
            return 373.15
        ln = __import__("math").log(self.pressure_pa / 101325.0)
        inv = 1.0 / 373.15 - (R * ln) / DH_VAP_WATER
        if inv <= 0:
            return 1e6
        return 1.0 / inv

    def water_phase(self) -> str:
        if not self.water:
            return "none (no water in this environment)"
        bp = self.water_boiling_point()
        if self.pressure_pa < 611.657 and self.temperature_k < 273.16:
            return "solid ice (below the triple point)"
        if self.temperature_k >= bp:
            return f"vapour (boils at {bp:.1f} K here)"
        if self.temperature_k <= 273.15:
            return "solid ice"
        return f"liquid (boils at {bp:.1f} K here)"

    def atmosphere_text(self) -> str:
        if not self.atmosphere:
            return "vacuum"
        parts = []
        for gas, frac in self.atmosphere:
            parts.append(f"{gas} {frac * 100:.2f}%")
        return " \u00b7 ".join(parts)

    def summary(self) -> str:
        return (
            f"{self.code} \u00b7 {self.name}\n"
            f"    temperature   {self.temperature_k:.2f} K ({self.temperature_c:.1f} \u00b0C)\n"
            f"    pressure      {self.pressure_pa / 1000:.3f} kPa "
            f"({self.pressure_atm:.3f} atm)\n"
            f"    atmosphere    {self.atmosphere_text()}\n"
            f"    water         {self.water_phase()}"
        )


# Reference conditions for each body: mean surface values, or the 1 bar level
# for the gas giants.  Atmosphere entries are mole fractions.
_BODY_TABLE: tuple[
    tuple[str, str, float, float, tuple[tuple[str, float], ...], bool, tuple[str, ...]], ...
] = (
    (
        "E1", "Earth surface (sea level, room temperature)", 288.15, 101325.0,
        (("N2", 0.7808), ("O2", 0.2095), ("Ar", 0.0093), ("CO2", 0.00042)), True,
        (
            "Open air at 1 atm; oxygen available for combustion.",
            "Liquid water available as solvent; gases can escape.",
        ),
    ),
    (
        "E2", "Earth orbit (outside the atmosphere)", 255.0, 0.0,
        (), False,
        (
            "Sunlit vacuum near Earth: about 390 K in sunlight, 200 K in shadow.",
            "No air, so no oxidation, no convection and no solvent.",
        ),
    ),
    (
        "S1", "Deep space (vacuum, cosmic background)", 2.725, 0.0,
        (), False,
        (
            "Hard vacuum: no oxidiser, no solvent, no convective cooling.",
            "Heat is lost only by radiation; any gas or vapour escapes.",
        ),
    ),
    (
        "ME1", "Mercury surface (day side)", 440.0, 1e-9,
        (), False,
        (
            "Essentially no atmosphere; 700 K by day, 100 K at night.",
            "No water, no oxidiser, no convection.",
        ),
    ),
    (
        "V1", "Venus surface", 737.0, 9.2e6,
        (("CO2", 0.965), ("N2", 0.035)), False,
        (
            "92 atm of carbon dioxide at 737 K: hot enough to melt lead.",
            "No liquid water, but CO2 is abundant and the pressure is enormous.",
        ),
    ),
    (
        "MO1", "Moon surface (average of day and night)", 250.0, 1e-9,
        (), False,
        (
            "Vacuum: 390 K by day, 100 K at night, no air and no water.",
        ),
    ),
    (
        "M1", "Mars surface", 210.0, 610.0,
        (("CO2", 0.951), ("N2", 0.028), ("Ar", 0.0197), ("O2", 0.0013)), False,
        (
            "Thin CO2 atmosphere at 0.006 atm; water is ice, not liquid.",
            "Dust storms and UV are the main environmental hazards.",
        ),
    ),
    (
        "J1", "Jupiter (1 bar level)", 165.0, 100000.0,
        (("H2", 0.898), ("He", 0.102)), False,
        (
            "Hydrogen atmosphere at 1 bar; hydrogen reduces most oxides.",
            "No water and no free oxygen.",
        ),
    ),
    (
        "SA1", "Saturn (1 bar level)", 134.0, 100000.0,
        (("H2", 0.963), ("He", 0.037)), False,
        (
            "Hydrogen atmosphere at 1 bar and 134 K: extremely cold.",
        ),
    ),
    (
        "T1", "Titan surface (moon of Saturn)", 94.0, 146700.0,
        (("N2", 0.95), ("CH4", 0.05)), False,
        (
            "Thick nitrogen atmosphere with methane, plus methane/ethane lakes.",
            "Water is rock-hard ice; the solvent here is methane.",
        ),
    ),
    (
        "U1", "Uranus (1 bar level)", 76.0, 100000.0,
        (("H2", 0.83), ("He", 0.15), ("CH4", 0.02)), False,
        (
            "Hydrogen atmosphere with methane: the coldest planetary atmosphere.",
        ),
    ),
    (
        "N1", "Neptune (1 bar level)", 72.0, 100000.0,
        (("H2", 0.80), ("He", 0.19), ("CH4", 0.01)), False,
        (
            "Hydrogen atmosphere with methane, extremely cold and windy.",
        ),
    ),
    (
        "P1", "Pluto surface", 44.0, 1.0,
        (("N2", 0.99), ("CH4", 0.01)), False,
        (
            "Very thin nitrogen atmosphere; water, methane and nitrogen are frozen.",
        ),
    ),
)

PRESETS = {
    code: Surroundings(
        code=code,
        name=name,
        temperature_k=temperature,
        pressure_pa=pressure,
        atmosphere=atmosphere,
        water=water,
        notes=notes,
    )
    for code, name, temperature, pressure, atmosphere, water, notes in _BODY_TABLE
}


def preset_table() -> list[tuple[str, str, float, float]]:
    """(code, name, temperature, pressure) for every built-in surroundings."""
    return [
        (env.code, env.name, env.temperature_k, env.pressure_pa)
        for env in PRESETS.values()
    ]


def get(code: str) -> Surroundings:
    key = code.strip().upper()
    if key in PRESETS:
        return PRESETS[key]
    raise KeyError(f"unknown surroundings {code!r}; try E1 or S1")


_PRESSURE_UNITS = {
    "pa": 1.0,
    "kpa": 1000.0,
    "mpa": 1e6,
    "atm": 101325.0,
    "bar": 100000.0,
    "mbar": 100.0,
    "torr": 101325.0 / 760.0,
    "mmhg": 101325.0 / 760.0,
    "psi": 6894.757,
}


def parse_surroundings(spec: str, base: Surroundings | None = None) -> Surroundings:
    """Parse ``E1``, ``S1`` or ``T=300K P=0.5atm O2=100% water=off``."""
    spec = spec.strip()
    if not spec:
        raise ValueError("empty surroundings specification")
    if spec.upper() in PRESETS:
        return PRESETS[spec.upper()]

    env = base or PRESETS["E1"]
    temperature = env.temperature_k
    pressure = env.pressure_pa
    atmosphere = dict(env.atmosphere)
    explicit_gases: dict[str, float] = {}
    water = env.water
    name = env.name

    for token in spec.replace(",", " ").split():
        if "=" not in token:
            raise ValueError(
                f"cannot read {token!r}; use key=value such as T=300K P=1atm O2=21%"
            )
        key, _, value = token.partition("=")
        key = key.strip().lower()
        value = value.strip()
        try:
            if key in ("t", "temp", "temperature"):
                if value.lower().endswith("c"):
                    temperature = float(value[:-1]) + 273.15
                elif value.lower().endswith("k"):
                    temperature = float(value[:-1])
                else:
                    temperature = float(value)
            elif key in ("p", "pressure"):
                lowered = value.lower()
                for unit, factor in sorted(_PRESSURE_UNITS.items(), key=lambda kv: -len(kv[0])):
                    if lowered.endswith(unit):
                        pressure = float(lowered[: -len(unit)]) * factor
                        break
                else:
                    pressure = float(lowered)
            elif key in ("water", "h2o"):
                water = value.lower() in ("on", "yes", "true", "1", "liquid", "available")
            elif key in ("name", "label"):
                name = value.replace("_", " ")
            elif key in ("co2", "n2", "o2", "ar", "he", "h2", "ch4", "nh3", "so2"):
                gas = key.upper()
                frac = float(value.rstrip("%")) / 100.0 if value.endswith("%") else float(value)
                if frac < 0:
                    raise ValueError("a gas fraction cannot be negative")
                explicit_gases[gas] = frac
                if frac <= 0:
                    atmosphere.pop(gas, None)
                else:
                    atmosphere[gas] = frac
            else:
                raise ValueError(f"unknown setting {key!r}")
        except ValueError as exc:
            raise ValueError(f"bad value for {key}: {value!r}") from exc

    if explicit_gases:
        # Gases the user named define the mix; the leftover background gases
        # share whatever room is left, so "O2=100%" means pure oxygen.
        chosen = sum(explicit_gases.values())
        if chosen > 1.0001:
            raise ValueError("atmosphere fractions add up to more than 100%")
        background = {
            gas: frac for gas, frac in atmosphere.items() if gas not in explicit_gases
        }
        background_total = sum(background.values())
        room = max(0.0, 1.0 - chosen)
        if background and background_total > 0 and room > 0:
            background = {
                gas: frac * room / background_total for gas, frac in background.items()
            }
        else:
            background = {}
        atmosphere = {**background, **{g: f for g, f in explicit_gases.items() if f > 0}}
    elif atmosphere:
        total = sum(atmosphere.values())
        if total > 1.0001:
            raise ValueError("atmosphere fractions add up to more than 100%")
        atmosphere = {gas: frac / total for gas, frac in atmosphere.items()} if total else {}
    if pressure > 0 and not atmosphere:
        atmosphere = {"N2": 0.7902, "O2": 0.2095, "Ar": 0.0003}

    return Surroundings(
        code="custom",
        name=name,
        temperature_k=temperature,
        pressure_pa=pressure,
        atmosphere=tuple(atmosphere.items()),
        water=water,
        notes=(f"custom surroundings derived from {env.code}",),
    )
