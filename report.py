"""Rendering an :class:`~maschem.solver.Analysis` for the terminal."""

from __future__ import annotations

import math
import shutil

from .solver import Analysis, time_series
from .thermo import cooling_temperature

WIDTH = 74
MIN_WIDTH = 46
MAX_WIDTH = 100


def _terminal_width() -> int:
    """Layout width that fits the terminal the report is printed to."""
    try:
        columns = shutil.get_terminal_size((80, 24)).columns
    except Exception:  # pragma: no cover - exotic terminals
        columns = 80
    return max(MIN_WIDTH, min(MAX_WIDTH, columns - 1))

RESET = "\x1b[0m"
BOLD = "\x1b[1m"
DIM = "\x1b[2m"
CYAN = "\x1b[36m"
YELLOW = "\x1b[33m"
GREEN = "\x1b[32m"
RED = "\x1b[31m"
MAGENTA = "\x1b[35m"


def num(value: float | None, sig: int = 4) -> str:
    if value is None:
        return "?"
    if value == 0:
        return "0"
    if abs(value) >= 1e5 or abs(value) < 1e-3:
        return f"{value:.{sig - 1}e}".replace("e+", "e").replace("e-0", "e-")
    return f"{value:.{sig}g}"


def _human_time(seconds: float) -> str:
    """A time that reads well on one line: 0.5 s, 45 s, 3 min, 2.4 h, 1.1 d."""
    try:
        value = float(seconds)
    except (TypeError, ValueError):
        return "?"
    if value < 1e-9:
        return "0 s"
    if value < 1.0:
        return f"{value:.2f} s"
    if value < 90.0:
        return f"{value:.0f} s"
    if value < 5400.0:
        return f"{value / 60.0:.0f} min"
    if value < 172800.0:
        return f"{value / 3600.0:.1f} h"
    return f"{value / 86400.0:.1f} d"


def _interpolate(points: list[tuple[float, float]], x: float) -> float | None:
    """Linear reading of a sampled curve, for drawing it between samples."""
    if not points:
        return None
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    for (xa, ya), (xb, yb) in zip(points, points[1:]):
        if xa <= x <= xb:
            if xb - xa < 1e-12:
                return yb
            return ya + (yb - ya) * ((x - xa) / (xb - xa))
    return points[-1][1]


def _plot(
    points: list[tuple[float, float | None]],
    height: int = 7,
    width: int | None = None,
) -> list[str]:
    """A small ASCII graph of a sampled series: value against time."""
    clean = [
        (float(x), float(y))
        for x, y in points
        if y is not None and math.isfinite(float(y))
    ]
    if len(clean) < 2:
        return []
    cells = width or max(16, min(WIDTH - 12, 52))
    first, last = clean[0][0], clean[-1][0]
    values = [y for _, y in clean]
    low, high = min(values), max(values)
    if high - low < 1e-9:
        low -= 0.5
        high += 0.5
    span = high - low
    decimals = 3 if span < 1 else 2 if span < 100 else 1 if span < 10000 else 0
    grid = [[" "] * cells for _ in range(height)]
    previous: int | None = None
    for column in range(cells):
        moment = first + (last - first) * column / max(cells - 1, 1)
        value = _interpolate(clean, moment)
        if value is None:
            previous = None
            continue
        row = int(round((high - value) / span * (height - 1)))
        if previous is not None:
            for between in range(min(previous, row), max(previous, row) + 1):
                if grid[between][column] == " ":
                    grid[between][column] = "\u2502"
        grid[row][column] = "\u00b7"
        previous = row
    lines: list[str] = []
    for index in range(height):
        value = high - span * index / (height - 1)
        label = (
            f"{value:>8.{decimals}f}"
            if index in (0, height // 2, height - 1)
            else " " * 8
        )
        lines.append((label + " \u2524" + "".join(grid[index])).rstrip())
    lines.append(" " * 8 + " \u2514" + "\u2500" * cells)
    half = max(cells // 2, 4)
    lines.append(" " * 10 + f"{_human_time(first):<{half}}{_human_time(last):>{half}}")
    return lines


def atoms_text(value: float | None) -> str:
    if value is None:
        return "?"
    return f"{value:.3e}".replace("e+", "e+")


class Report:
    def __init__(self, colour: bool = True):
        self.lines: list[str] = []
        self.colour = colour

    def paint(self, text: str, code: str) -> str:
        return f"{code}{text}{RESET}" if self.colour else text

    def blank(self) -> None:
        self.lines.append("")

    def banner(self, env) -> None:
        self.lines.append(self.paint("MAS CHEM", BOLD + CYAN) + self.paint("  \u00b7  terminal reaction solver", DIM))
        label = f"surroundings {env.code} \u00b7 {env.name}"
        for chunk in _wrap(label, WIDTH):
            self.lines.append(self.paint(chunk, DIM))

    def section(self, title: str) -> None:
        label = f" {title} "
        pad = max(WIDTH - len(label) - 2, 0)
        self.lines.append(self.paint("\u2500" + label + "\u2500" * pad, CYAN))

    def row(self, label: str, value: str, width: int = 14) -> None:
        prefix = f"  {label:<{width}}"
        text = prefix + value
        if len(text) <= WIDTH:
            self.lines.append(text)
            return
        # Long values (equations, atmospheres) wrap onto aligned continuations
        # instead of being cut by the terminal.
        available = max(12, WIDTH - len(prefix))
        for index, chunk in enumerate(_wrap(value, available)):
            self.lines.append(prefix + chunk if index == 0 else " " * len(prefix) + chunk)

    def series_table(self, rows: list[dict]) -> None:
        """Time, temperature and pH on the same line for every sample."""
        has_ph = any(row.get("ph") is not None for row in rows)
        head = f"  {'t/s':>8}  {'T/K':>7}  {'\u00b0C':>7}"
        if has_ph:
            head += f"  {'pH':>6}"
        self.lines.append(self.paint(head, BOLD))
        for row in rows:
            line = f"  {row['t']:>8.1f}  {row['T']:>7.1f}  {row['C']:>7.1f}"
            if has_ph:
                ph = row.get("ph")
                line += f"  {ph:>6.2f}" if ph is not None else f"  {'-':>6}"
            self.lines.append(line)

    def note(self, text: str, mark: str = "\u00b7") -> None:
        for index, chunk in enumerate(_wrap(text, WIDTH - 8)):
            prefix = f"    {mark} " if index == 0 else "      "
            self.lines.append(self.paint(prefix + chunk, DIM))

    def render(self, analysis: Analysis) -> str:
        env = analysis.env
        self.banner(env)

        if analysis.error:
            self.blank()
            self.section("NO RESULT")
            self.lines.append(self.paint(f"  {analysis.error}", YELLOW))
            return "\n".join(self.lines)

        self.blank()
        self.section("SURROUNDINGS")
        self.row("temperature", f"{env.temperature_k:.2f} K  ({env.temperature_c:.1f} \u00b0C)")
        self.row("pressure", f"{env.pressure_pa / 1000:.3f} kPa  ({env.pressure_atm:.3f} atm)")
        self.row("atmosphere", env.atmosphere_text(), width=14)
        self.row("water", env.water_phase())

        self.blank()
        self.section("REACTION")
        self.row("input", analysis.input_text.strip())
        self.row("source", analysis.source)
        self.row("balanced", analysis.equation)
        if analysis.limiting:
            self.row("limiting", analysis.limiting)

        self.blank()
        self.section("ENERGY & TEMPERATURE")
        if analysis.dh_per is not None:
            kind = "exothermic" if analysis.dh_per < 0 else "endothermic"
            self.row("enthalpy", f"{analysis.dh_per:+.1f} kJ per mol of reaction  ({kind})")
            if analysis.dh_total is not None and not analysis.per_mole:
                direction = "released" if analysis.dh_total < 0 else "absorbed"
                self.row("heat", f"{abs(analysis.dh_total):.2f} kJ {direction}")
        else:
            self.row("enthalpy", "not available (missing formation data)")

        self.row("t = 0", f"{analysis.t_initial:.2f} K  (surroundings)")
        if analysis.t_peak is not None:
            delta = analysis.t_peak - analysis.t_initial
            if delta < -0.5:
                self.row(
                    "dip",
                    f"{analysis.t_peak:.0f} K  (adiabatic estimate, {delta:+.0f} K; "
                    "the reaction absorbs heat)",
                )
            else:
                self.row(
                    "peak",
                    f"{analysis.t_peak:.0f} K  (adiabatic estimate, {delta:+.0f} K)",
                )
        else:
            self.row("peak", "not available")
        self.row("settled", f"{analysis.t_final:.2f} K  (after heat exchange with surroundings)")
        if analysis.peak_phases:
            self.row("phases", analysis.peak_phases[0])
            for extra in analysis.peak_phases[1:]:
                self.row("", extra)

        if analysis.timeseries:
            self.blank()
            self.section("TEMPERATURE & pH VS TIME")
            loss = analysis.loss_coefficient or 0.5
            self.row(
                "model",
                f"Newton cooling, \u03c4 = {analysis.cooling_tau:.3g} s "
                f"({analysis.heat_capacity:.3g} J/K, {loss:g} W/K loss to the room)",
            )
            if analysis.curve_window:
                self.row(
                    "window",
                    f"{_human_time(analysis.curve_window)} (4\u03c4), sampled every "
                    f"{_human_time(analysis.sample_step or 0.0)}",
                )
            self.series_table(analysis.timeseries)
            if analysis.temperature_curve:
                self.lines.append(self.paint("  temperature", BOLD))
                self.lines.extend(_plot(analysis.temperature_curve))
            if len(analysis.ph_curve) >= 2:
                self.lines.append(self.paint("  pH", BOLD))
                self.lines.extend(_plot(analysis.ph_curve))
                self.note(
                    "the mixture is fully reacted at t = 0; pH moves because Kw "
                    "and the acid equilibria shift with temperature"
                )

        if analysis.initial_ph is not None or analysis.final_ph is not None:
            self.blank()
            self.section("pH")
            if analysis.initial_ph is not None and analysis.initial_ph.ph is not None:
                self.row("before", f"{analysis.initial_ph.ph:.2f}")
            else:
                self.row("before", "no aqueous reactants")
            if analysis.final_ph is not None and analysis.final_ph.ph is not None:
                self.row("after", f"{analysis.final_ph.ph:.2f}")
            else:
                self.row("after", "no aqueous products")
            neutral = (analysis.final_ph or analysis.initial_ph).neutral_ph
            self.note(
                f"neutral pH at {analysis.t_initial:.1f} K is {neutral:.2f}; "
                f"assumes {analysis.volume_l:g} L of water"
            )
            for detail in (analysis.final_ph.details if analysis.final_ph else []):
                self.note(detail)
        elif env.water:
            self.blank()
            self.section("pH")
            self.row("result", "no dissolved acid/base present")

        self.blank()
        self.section("COLOUR")
        self.row("before", analysis.colour_before or "colourless solution")
        self.row("after", _product_colours(analysis))

        self.blank()
        self.section("REDOX")
        red = analysis.redox
        if red is None:
            self.row("result", "not analysed")
        elif not red.happened:
            self.row("result", "no \u2014 no element changes oxidation state")
        else:
            self.row("result", self.paint("yes", GREEN))
            for change in red.changes:
                verb = "oxidised" if change.oxidised else "reduced"
                arrow = f"{_ox(change.before)} \u2192 {_ox(change.after)}"
                sign = "loses" if change.oxidised else "gains"
                detail = f"{verb}, {change.atoms:g} atoms"
                value = (
                    f"{arrow:<12}{sign} {abs(change.per_atom):.3g} e- per atom  ({detail})"
                )
                if change.transition:
                    value += f"  {change.transition}"
                self.row(f"  {change.element}", value)
            total_electrons = red.electrons * analysis.extent
            self.row("electrons", f"{num(red.electrons, 4)} per reaction")
            self.row(
                "  total",
                f"{num(total_electrons, 4)} mol = {atoms_text(total_electrons * 6.02214076e23)} electrons",
            )
            oxidised = ", ".join(red.oxidised_elements()) or "\u2014"
            reduced = ", ".join(red.reduced_elements()) or "\u2014"
            self.row("reducing", oxidised)
            self.row("oxidising", reduced)

        self.blank()
        self.section("PRODUCTS")
        # Narrow terminals lose the atom-count column rather than wrapping it.
        show_atoms = WIDTH >= 66
        name_width = 20 if show_atoms else 14
        header = f"  {'species':<{name_width}}{'phase':<7}{'mol':>8}{'mass/g':>10}"
        if show_atoms:
            header += f"{'atoms':>14}"
        self.lines.append(self.paint(header, BOLD))
        for row in analysis.product_rows:
            sp = row["species"]
            cells = (
                f"  {sp.display:<{name_width}}{sp.state_text:<7}"
                f"{num(row['moles']):>8}{num(row['mass'], 4):>10}"
            )
            if show_atoms:
                cells += f"{atoms_text(row['atoms']):>14}"
            self.lines.append(cells)
        for row in analysis.unreacted:
            sp = row["species"]
            cells = (
                f"  {sp.display + ' (left over)':<{name_width}}{sp.state_text:<7}"
                f"{num(row['moles']):>8}{num(row['mass'], 4):>10}"
            )
            self.lines.append(self.paint(cells, DIM))
        self.row("mass balance", f"{num(analysis.mass_in, 4)} g in \u2192 {num(analysis.mass_out, 4)} g out")

        self.blank()
        self.section("ATOMS")
        if analysis.element_totals:
            for symbol in sorted(analysis.element_totals):
                entry = analysis.element_totals[symbol]
                atoms = atoms_text(entry["atoms"])
                moles = num(entry["moles"], 4)
                mass = num(entry["mass"], 4)
                if WIDTH >= 66:
                    self.lines.append(
                        f"  {symbol:<4}{atoms:>14} atoms{moles:>12} mol{mass:>10} g"
                    )
                else:
                    self.lines.append(
                        f"  {symbol:<4}{atoms:>12} {moles:>8} mol {mass:>8} g"
                    )
        else:
            self.row("result", "no atoms counted")

        notes = list(analysis.notes)
        if analysis.warnings:
            notes += analysis.warnings
        if notes:
            self.blank()
            self.section("NOTES")
            for text in notes:
                self.note(text)
        return "\n".join(self.lines)


def _wrap(text: str, width: int) -> list[str]:
    words = text.split()
    if not words:
        return [""]
    lines: list[str] = []
    current = words[0]
    for word in words[1:]:
        if len(current) + 1 + len(word) <= width:
            current += " " + word
        else:
            lines.append(current)
            current = word
    lines.append(current)
    return lines


def _ox(value: float) -> str:
    if abs(value) < 0.005:
        return "0"
    if abs(value - round(value)) < 0.01:
        value = round(value)
        return f"{value:+d}"
    return f"{value:+.2f}"


def _product_colours(analysis: Analysis) -> str:
    from .data import COMPOUND_COLOURS, GAS_COLOURS

    if analysis.final_ph is not None:
        from .solver import _solution_colour

        aqueous = []
        for row in analysis.product_rows:
            sp = row["species"]
            if sp.state == "aq":
                aqueous.append((sp, row["moles"]))
        for row in analysis.unreacted:
            if row["species"].state == "aq":
                aqueous.append((row["species"], row["moles"]))
        solution = _solution_colour(aqueous)
        solids = [
            COMPOUND_COLOURS.get(row["species"].formula)
            for row in analysis.product_rows
            if row["species"].state in ("s", "l")
            and COMPOUND_COLOURS.get(row["species"].formula)
        ]
        parts = []
        if solution != "colourless solution":
            parts.append(solution)
        if solids:
            parts.append("; ".join(solids))
        if parts:
            return " \u00b7 ".join(parts)
        if aqueous:
            return "colourless solution"

    parts = []
    for row in analysis.product_rows:
        sp = row["species"]
        if sp.state == "aq":
            continue
        if sp.state == "g" and sp.formula in GAS_COLOURS:
            parts.append(f"{sp.display} is {GAS_COLOURS[sp.formula]}")
        elif sp.formula in COMPOUND_COLOURS:
            parts.append(f"{sp.display} is {COMPOUND_COLOURS[sp.formula]}")
    return "; ".join(parts) if parts else "no colour data"


def curve_table(
    analysis: Analysis, step: float | None = None, end: float | None = None
) -> str:
    """Temperature and pH at a chosen resolution, for the ``curve`` command."""
    global WIDTH
    WIDTH = _terminal_width()
    if analysis.t_peak is None or not analysis.cooling_tau:
        return (
            "no temperature curve: this reaction has no available enthalpy "
            "(or it is endothermic, so it never heats up)"
        )
    rows = time_series(analysis, step, end)
    if not rows:
        return "no temperature curve for this reaction"
    tau = analysis.cooling_tau
    used = rows[1]["t"] - rows[0]["t"] if len(rows) > 1 else 0.0
    toward = "dip" if analysis.t_peak < analysis.t_initial else "peak"
    lines = [
        f"temperature and pH every {_human_time(used)} for {_human_time(rows[-1]['t'])} "
        f"(\u03c4 = {tau:.3g} s, {toward} {analysis.t_peak:.0f} K)"
    ]
    table = Report(colour=False)
    table.series_table(rows)
    lines.extend(table.lines)
    lines.append("  temperature")
    lines.extend(_plot([(row["t"], row["T"]) for row in rows]))
    ph_points = [(row["t"], row["ph"]) for row in rows if row["ph"] is not None]
    if len(ph_points) >= 2:
        lines.append("  pH")
        lines.extend(_plot(ph_points))
    return "\n".join(lines)


def render(analysis: Analysis, colour: bool = True) -> str:
    global WIDTH
    WIDTH = _terminal_width()
    return Report(colour=colour).render(analysis)
