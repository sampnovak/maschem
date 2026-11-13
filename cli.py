"""Interactive terminal front-end for MasChem."""

from __future__ import annotations

import os
import re
import sys

from .data import DHF
from .formula import FormulaError, parse_species
from .report import BOLD, CYAN, DIM, GREEN, RESET, YELLOW, curve_table, render
from .solver import Analysis, analyse
from .surroundings import (
    PRESETS,
    Surroundings,
    get,
    parse_surroundings,
    preset_table,
)
from .thermo import LOSS_COEFFICIENT

EXAMPLES = [
    "2g Fe + 3g S",
    "2g of iron mixed with 2g of sulfur",
    "CH4 + 2O2 -> CO2 + 2H2O",
    "25g CaCO3 -> CaO + CO2",
    "Fe + CuSO4 -> FeSO4 + Cu",
    "100g Al + 100g Fe2O3 @S1",
    "0.5 mol HCl + 0.5 mol NaOH",
    "10g Mg + 10g O2",
]

HELP = """\
MasChem reads either a balanced-style equation or a list of amounts.

  equation     CH4 + 2O2 -> CO2 + 2H2O
  amounts      2g Fe + 3g S
              2g of iron mixed with 2g of sulfur
              0.5 mol HCl + 0.5 mol NaOH
              10g Mg + 10g O2

If you leave out the products, MasChem predicts them from its reaction
database and inorganic rules (combustion, neutralisation, precipitation,
displacement, synthesis and decomposition).

Amounts accept g, mg, kg, mol and mmol.  Without amounts the answer is
per mole of the balanced reaction.

Commands
  surroundings              show the current environment
  surroundings list         every built-in environment
  surroundings E1           Earth surface (288 K, 1 atm air, water)
  surroundings E2           Earth orbit (255 K, vacuum)     S1 deep space
  surroundings M1           Mars (210 K, 0.006 atm CO2)
  surroundings V1           Venus   J1 Jupiter   SA1 Saturn   T1 Titan
  surroundings ME1          Mercury MO1 Moon      U1 Uranus    N1 Neptune
  surroundings P1           Pluto
  surroundings T=350K P=0.5atm O2=100% water=off
  volume 2                  solution volume: sets pH strength AND thermal mass
                            (2 L heats and cools twice as slowly as 1 L)
  cooling 0.5               heat lost to the room, W/K (0.05 dewar, 5 fan)
  dhf MnCl2(aq) -555.1      add a formation enthalpy for your own compound
  curve                     temperature and pH over four time constants
  curve 1 120               ... every second for 120 s (last reaction)
  examples                  show example inputs
  help                      this text
  quit / exit               leave

Append @E1, @M1, @S1 ... to one line to run it in that environment, e.g.
  10g Mg + 10g O2 @S1

Run the built-in solver checks with
  python -m maschem --selftest

Full guide: README.md, docs/usage.md and docs/how-it-works.md in the maschem folder.
"""


def _enable_ansi() -> None:
    if os.name != "nt":
        return
    try:
        import ctypes

        kernel32 = ctypes.windll.kernel32
        handle = kernel32.GetStdHandle(-11)
        mode = ctypes.c_uint32()
        if kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            kernel32.SetConsoleMode(handle, mode.value | 0x0004)
    except Exception:
        pass


def _supports_colour() -> bool:
    if os.environ.get("NO_COLOR"):
        return False
    return sys.stdout.isatty()


def _resolve_env(spec: str, base: Surroundings | None = None) -> Surroundings:
    spec = spec.strip()
    if not spec:
        raise ValueError("empty surroundings specification")
    if any(mark in spec for mark in "=,"):
        return parse_surroundings(spec, base)
    try:
        return get(spec)
    except KeyError:
        raise ValueError(f"unknown surroundings {spec!r}; try E1 or S1") from None


def _extract_surroundings(
    text: str, base: Surroundings | None = None
) -> tuple[str, Surroundings | None]:
    match = re.search(r"@\s*([^@]*?)\s*$", text)
    if not match or not match.group(1).strip():
        return text, None
    return text[: match.start()].strip(), _resolve_env(match.group(1), base)


class Shell:
    def __init__(self) -> None:
        self.env = PRESETS["E1"]
        self.volume = 1.0
        self.colour = _supports_colour()
        self.loss = LOSS_COEFFICIENT
        self.last: Analysis | None = None

    def run_once(self, text: str) -> str:
        body, override = _extract_surroundings(text, self.env)
        env = override or self.env
        analysis = analyse(body, env, self.volume, loss=self.loss)
        return render(analysis, colour=self.colour)

    def repl(self) -> None:
        _enable_ansi()
        print(self._paint("MAS CHEM", BOLD + CYAN) + self._paint("  \u00b7  terminal reaction solver", DIM))
        print(self._paint("surroundings E1 \u00b7 Earth surface (15 \u00b0C, 1 atm air, liquid water)", DIM))
        print(self._paint("type `help` for commands, or just give me a reaction", DIM))
        while True:
            try:
                line = input(self._paint("\nmaschem> ", GREEN)).strip()
            except (EOFError, KeyboardInterrupt):
                print()
                break
            if not line:
                continue
            command = line.split()[0].lower()
            if command in ("quit", "exit", "q"):
                break
            if command in ("help", "?", "h"):
                print(HELP)
                continue
            if command == "examples":
                for example in EXAMPLES:
                    print(f"  {example}")
                continue
            if command in ("surroundings", "env"):
                rest = line[len(command) :].strip()
                if not rest:
                    print(self.env.summary())
                elif rest.lower() in ("list", "all", "bodies"):
                    print(f"  {'code':<4} {'T/K':>7}  {'pressure':>10}  name")
                    for code, name, temp, pressure in preset_table():
                        bar = "vacuum" if pressure < 1.0 else f"{pressure / 101325:.3g} atm"
                        print(f"  {code:<4} {temp:>7.2f}  {bar:>10}  {name}")
                else:
                    try:
                        self.env = parse_surroundings(rest, self.env)
                        print(self.env.summary())
                    except (ValueError, KeyError) as exc:
                        print(self._paint(f"! {exc}", YELLOW))
                continue
            if command in ("volume", "vol"):
                rest = line[len(command) :].strip().lower().rstrip("l")
                try:
                    value = float(rest)
                    if value <= 0:
                        raise ValueError
                    self.volume = value
                    print(f"solution volume set to {value:g} L")
                except ValueError:
                    print(self._paint("! usage: volume 1.5", YELLOW))
                continue
            if command in ("cooling", "loss"):
                rest = line[len(command) :].strip()
                if not rest:
                    print(
                        f"heat loss to the room: {self.loss:g} W/K "
                        "(0.05 dewar, 0.5 open beaker, 5 fan-cooled)"
                    )
                    continue
                try:
                    value = float(rest.split()[0])
                    if value <= 0:
                        raise ValueError
                except ValueError:
                    print(self._paint("! usage: cooling 0.5   (watts per kelvin)", YELLOW))
                    continue
                self.loss = value
                print(f"heat loss set to {value:g} W/K")
                continue
            if command in ("dhf", "enthalpy"):
                rest = line[len(command) :].split()
                if not rest:
                    print(self._paint("! usage: dhf MnCl2(aq) -555.1", YELLOW))
                    continue
                try:
                    species = parse_species(rest[0])
                except FormulaError as exc:
                    print(self._paint(f"! {exc}", YELLOW))
                    continue
                key = (
                    f"{species.formula}({species.state})"
                    if species.state
                    else species.formula
                )
                if len(rest) == 1:
                    known = DHF.get(key, DHF.get(species.formula))
                    if known is None:
                        print(self._paint(f"! no formation enthalpy for {key}", YELLOW))
                    else:
                        print(f"{key}: {known:+.1f} kJ/mol")
                    continue
                try:
                    value = float(rest[1])
                except ValueError:
                    print(self._paint("! the enthalpy must be a number in kJ/mol", YELLOW))
                    continue
                DHF[key] = value
                print(f"{key}: {value:+.1f} kJ/mol  (this session only)")
                continue
            if command in ("curve", "temp"):
                rest = line[len(command) :].split()
                try:
                    step = float(rest[0]) if rest else None
                    end = float(rest[1]) if len(rest) > 1 else None
                    if (step is not None and step <= 0) or (
                        end is not None and end <= 0
                    ):
                        raise ValueError
                except ValueError:
                    print(self._paint("! usage: curve [step-seconds] [end-seconds]", YELLOW))
                    continue
                if self.last is None:
                    print(self._paint("! run a reaction first, then `curve` shows its temperature", YELLOW))
                else:
                    print(curve_table(self.last, step, end))
                continue
            if command in ("clear", "cls"):
                os.system("cls" if os.name == "nt" else "clear")
                continue

            try:
                body, override = _extract_surroundings(line, self.env)
                env = override or self.env
                analysis = analyse(body, env, self.volume, loss=self.loss)
                self.last = analysis
                print()
                print(render(analysis, colour=self.colour))
            except Exception as exc:  # keep the REPL alive
                print(self._paint(f"! {exc}", YELLOW))

    def _paint(self, text: str, code: str) -> str:
        return f"{code}{text}{RESET}" if self.colour else text


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    shell = Shell()
    env_override = None
    expression: list[str] = []

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg in ("-e", "--env", "--surroundings") and i + 1 < len(argv):
            try:
                env_override = _resolve_env(argv[i + 1])
            except ValueError as exc:
                print(exc)
                return 2
            i += 2
            continue
        if arg in ("--cooling", "--loss") and i + 1 < len(argv):
            try:
                value = float(argv[i + 1])
                if value <= 0:
                    raise ValueError
            except ValueError:
                print("cooling must be a positive number of W/K")
                return 2
            shell.loss = value
            i += 2
            continue
        if arg in ("-v", "--volume") and i + 1 < len(argv):
            try:
                shell.volume = float(argv[i + 1])
            except ValueError:
                print("volume must be a number")
                return 2
            i += 2
            continue
        if arg in ("-h", "--help"):
            print(HELP)
            return 0
        if arg == "--no-colour":
            shell.colour = False
            i += 1
            continue
        if arg in ("--selftest", "--self-test"):
            from .selftest import main as run_checks

            return run_checks()
        expression.append(arg)
        i += 1

    if expression:
        text = " ".join(expression)
        if text.upper() in PRESETS:
            print(PRESETS[text.upper()].summary())
            return 0
        try:
            body, inline = _extract_surroundings(text, env_override or shell.env)
        except ValueError as exc:
            print(exc)
            return 2
        env = inline or env_override or shell.env
        analysis = analyse(body, env, shell.volume, loss=shell.loss)
        print(render(analysis, colour=shell.colour))
        return 0 if analysis.error is None else 1

    if env_override is not None:
        shell.env = env_override
    shell.repl()
    return 0
