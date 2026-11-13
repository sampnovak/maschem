"""Built-in checks for MasChem: ``python -m maschem --selftest``.

Each case states the balanced equation it expects, so a regression in the
reaction database, the heuristics, the balancer or the report shows up
immediately.  Invariants (mass balance, sane amounts, sensible pH) are checked
for every successful analysis, and a few checks cover the redox bookkeeping
and the temperature curve.
"""

from __future__ import annotations

import sys
import time

from .report import curve_table, render
from .solver import analyse
from .surroundings import PRESETS

ENV = PRESETS["E1"]

# (input, expected balanced equation, expected prediction source)
CASES: list[tuple[str, str, str]] = [
    # --- reaction database (must stay first) -----------------------------
    ("2g Fe + 3g S", "Fe(s) + S(s) → FeS(s)", "database"),
    ("CH4 + 2O2 -> CO2 + 2H2O", "CH4(g) + 2 O2(g) → CO2(g) + 2 H2O(l)", "given"),
    ("Na + H2O", "2 Na(s) + 2 H2O(l) → 2 NaOH(aq) + H2(g)", "database"),
    ("Fe + CuSO4 -> FeSO4 + Cu", "Fe(s) + CuSO4(aq) → FeSO4(aq) + Cu(s)", "given"),
    ("vinegar + baking soda",
     "CH3COOH(aq) + NaHCO3(aq) → CH3COONa(aq) + H2O(l) + CO2(g)", "database"),
    ("10g Mg + 10g O2", "2 Mg(s) + O2(g) → 2 MgO(s)", "database"),
    ("25g CaCO3 -> CaO + CO2", "CaCO3(s) → CaO(s) + CO2(g)", "given"),
    ("100g Al + 100g Fe2O3", "2 Al(s) + Fe2O3(s) → Al2O3(s) + 2 Fe(s)", "database"),
    # --- metal oxide + acid ---------------------------------------------
    ("Fe2O3 + HCl", "Fe2O3(s) + 6 HCl(aq) → 2 FeCl3(aq) + 3 H2O(l)",
     "metal oxide + acid"),
    ("CuO + H2SO4", "CuO(s) + H2SO4(aq) → CuSO4(aq) + H2O(l)",
     "metal oxide + acid"),
    ("ZnO + HNO3", "ZnO(s) + 2 HNO3(aq) → Zn(NO3)2(aq) + H2O(l)",
     "metal oxide + acid"),
    ("Al2O3 + HCl", "Al2O3(s) + 6 HCl(aq) → 2 AlCl3(aq) + 3 H2O(l)",
     "metal oxide + acid"),
    ("CaO + H3PO4", "3 CaO(s) + 2 H3PO4(aq) → Ca3(PO4)2(s) + 3 H2O(l)",
     "metal oxide + acid"),
    # --- acidified oxidiser ---------------------------------------------
    ("KMnO4 + HCl",
     "2 KMnO4(aq) + 16 HCl(aq) → 2 KCl(aq) + 2 MnCl2(aq) + 8 H2O(l) + 5 Cl2(g)",
     "oxidiser"),
    ("K2Cr2O7 + HCl",
     "K2Cr2O7(aq) + 14 HCl(aq) → 2 KCl(aq) + 2 CrCl3(aq) + 7 H2O(l) + 3 Cl2(g)",
     "oxidiser"),
    ("MnO2 + HCl", "MnO2(s) + 4 HCl(aq) → MnCl2(aq) + 2 H2O(l) + Cl2(g)",
     "oxidiser"),
    ("PbO2 + HCl", "PbO2(s) + 4 HCl(aq) → PbCl2(s) + 2 H2O(l) + Cl2(g)",
     "oxidiser"),
    ("KMnO4 + HBr",
     "2 KMnO4(aq) + 16 HBr(aq) → 2 KBr(aq) + 2 MnBr2(aq) + 8 H2O(l) + 5 Br2(l)",
     "oxidiser"),
    # --- acidic oxide + base --------------------------------------------
    ("CO2 + NaOH", "CO2(g) + 2 NaOH(aq) → Na2CO3(aq) + H2O(l)",
     "acidic oxide + base"),
    ("SO2 + NaOH", "SO2(g) + 2 NaOH(aq) → Na2SO3(aq) + H2O(l)",
     "acidic oxide + base"),
    ("SO3 + KOH", "SO3(s) + 2 KOH(aq) → K2SO4(aq) + H2O(l)",
     "acidic oxide + base"),
    ("P4O10 + NaOH", "P4O10(s) + 12 NaOH(aq) → 4 Na3PO4(aq) + 6 H2O(l)",
     "acidic oxide + base"),
    # --- acid + sulfite / sulfide ---------------------------------------
    ("Na2SO3 + HCl", "Na2SO3(aq) + 2 HCl(aq) → 2 NaCl(aq) + H2O(l) + SO2(g)",
     "sulfite"),
    ("FeS + HCl", "FeS(s) + 2 HCl(aq) → FeCl2(aq) + H2S(g)", "sulfide"),
    ("Na2S + H2SO4", "Na2S(aq) + H2SO4(aq) → Na2SO4(aq) + H2S(g)", "sulfide"),
    # --- oxide + hydrogen -----------------------------------------------
    ("Fe2O3 + H2", "Fe2O3(s) + 3 H2(g) → 2 Fe(s) + 3 H2O(l)", "hydrogen"),
    ("NiO + H2", "NiO(s) + H2(g) → Ni(s) + H2O(l)", "hydrogen"),
    # --- single-reactant decompositions ---------------------------------
    ("MgCO3", "MgCO3(s) → MgO(s) + CO2(g)", "carbonate decomposition"),
    ("ZnCO3", "ZnCO3(s) → ZnO(s) + CO2(g)", "carbonate decomposition"),
    ("Cu(OH)2", "Cu(OH)2(s) → CuO(s) + H2O(l)", "hydroxide decomposition"),
    ("Fe(OH)3", "2 Fe(OH)3(s) → Fe2O3(s) + 3 H2O(l)", "hydroxide decomposition"),
    # --- acids, bases, pH and redox in solution --------------------------
    ("NH4+ + OH- -> NH3 + H2O", "NH4+(aq) + OH−(aq) → NH3(g) + H2O(l)", "given"),
    ("1 mol Na2CO3 + 1 mol CH3COOH",
     "Na2CO3(aq) + 2 CH3COOH(aq) → 2 NaCH3COO(aq) + H2O(l) + CO2(g)",
     "acid + carbonate"),
    ("1 mol FeCl3 + 3 mol NaOH",
     "FeCl3(aq) + 3 NaOH(aq) → Fe(OH)3(s) + 3 NaCl(aq)", "database"),
    ("1 mol H2SO4 + 2 mol NaOH",
     "H2SO4(aq) + 2 NaOH(aq) → Na2SO4(aq) + 2 H2O(l)", "neutralisation"),
    ("1 mol KOH + 1 mol HNO3",
     "KOH(aq) + HNO3(aq) → KNO3(aq) + H2O(l)", "neutralisation"),
    ("Zn + CuSO4", "Zn(s) + CuSO4(aq) → ZnSO4(aq) + Cu(s)", "database"),
    ("10g NaCl + 10g AgNO3",
     "NaCl(aq) + AgNO3(aq) → AgCl(s) + NaNO3(aq)", "database"),
    ("Cl2 + NaBr", "Cl2(g) + 2 NaBr(aq) → 2 NaCl(aq) + Br2(l)", "database"),
    # --- rule-based, not a fixed list: any carbonate salt works ----------
    ("Li2CO3 + HCl", "Li2CO3(aq) + 2 HCl(aq) → 2 LiCl(aq) + H2O(l) + CO2(g)",
     "acid + carbonate"),
]

# Inputs that must fail cleanly rather than invent products.
ERRORS = [
    "Cu + HCl",                      # below hydrogen in the activity series
    "garbage text",                  # unparseable
    "",                              # empty
    "-> FeS",                        # no reactants
    "AgNO3 + NaCl -> AgNO3",         # a reactant would be left out
    "Al2O3 + NaOH",                  # amphoteric: not modelled, must not guess
    "ZnO + H2",                      # hydrogen cannot reduce zinc oxide
    "CaCO3 + Fe2O3",                 # no sensible reaction
]


def _curve_rows(text: str) -> int:
    """Data rows of a temperature/pH table (not headers, not graph axes)."""
    count = 0
    for line in text.splitlines():
        fields = line.split()
        if "\u2524" in line or len(fields) < 3:
            continue
        try:
            float(fields[0])
            float(fields[1])
        except ValueError:
            continue
        count += 1
    return count


def _model_checks() -> list[tuple[str, bool]]:
    permanganate = analyse("KMnO4 + HCl", ENV, 1.0)
    fast = analyse("10g Mg + 10g O2", ENV, 1.0, loss=5.0)
    slow = analyse("10g Mg + 10g O2", ENV, 1.0, loss=0.05)
    report_text = render(permanganate, colour=False)
    return [
        ("curve window scales with the time constant",
         bool(fast.curve_window) and bool(slow.curve_window)
         and fast.curve_window * 50 < slow.curve_window),
        ("permanganate window spans hours, not seconds",
         (permanganate.curve_window or 0.0) > 3600.0
         and len(permanganate.timeseries) > 12),
        ("window really shows the dip recovering",
         bool(permanganate.timeseries)
         and abs(permanganate.timeseries[-1]["T"] - ENV.temperature_k) < 1.0),
        ("pH reported on the same timetable as temperature",
         len(permanganate.ph_curve) > 5
         and abs(permanganate.ph_curve[0][1] - permanganate.ph_curve[-1][1]) > 0.05),
        ("report draws a temperature and a pH graph", report_text.count("\u2524") >= 2),
        ("permanganate has a formation-enthalpy temperature",
         permanganate.dh_per is not None and permanganate.t_peak is not None),
        ("heat-absorbing reaction reported as a dip",
         permanganate.t_peak is not None
         and permanganate.t_peak < ENV.temperature_k),
        ("permanganate gets a temperature curve",
         len(permanganate.temperature_curve) > 5),
        ("cooling coefficient changes the time constant",
         bool(fast.cooling_tau) and bool(slow.cooling_tau)
         and fast.cooling_tau < slow.cooling_tau),
        ("cooling coefficient is recorded on the analysis",
         abs((fast.loss_coefficient or 0.0) - 5.0) < 1e-9),
        ("solar-system presets are present", len(PRESETS) >= 13),
        ("Mars is thin and cold",
         "M1" in PRESETS
         and PRESETS["M1"].temperature_k < 220
         and PRESETS["M1"].pressure_atm < 0.01),
        ("Earth orbit is a vacuum", "E2" in PRESETS and PRESETS["E2"].is_vacuum),
        ("Venus is hot and thick",
         PRESETS["V1"].temperature_k > 700
         and PRESETS["V1"].pressure_atm > 90),
    ]


def _redox_checks() -> list[tuple[str, bool]]:
    result = analyse("KMnO4 + HCl", ENV, 1.0)
    changes = result.redox.changes if result.redox else []
    chlorine = [c for c in changes if c.element == "Cl"]
    manganese = [c for c in changes if c.element == "Mn"]
    return [
        ("chlorine reported per species as -1 -> 0",
         len(chlorine) == 1 and chlorine[0].before == -1 and chlorine[0].after == 0),
        ("all ten chloride atoms accounted for",
         len(chlorine) == 1 and abs(chlorine[0].atoms - 10) < 1e-9),
        ("chlorine species named (HCl -> Cl2)",
         len(chlorine) == 1
         and chlorine[0].from_species == "HCl"
         and chlorine[0].to_species == "Cl2"),
        ("manganese reduced +7 -> +2",
         len(manganese) == 1 and manganese[0].before == 7 and manganese[0].after == 2),
        ("ten electrons transferred",
         result.redox is not None and abs(result.redox.electrons - 10) < 1e-9),
    ]


def _temperature_checks() -> list[tuple[str, bool]]:
    result = analyse("10g Mg + 10g O2", ENV, 1.0)
    curve = result.temperature_curve
    peak = result.t_peak or 0.0
    return [
        ("temperature curve is sampled", len(curve) > 5),
        ("curve starts at the adiabatic peak at t = 0",
         bool(curve) and abs(curve[0][0]) < 1e-9 and abs(curve[0][1] - peak) < 1e-6),
        ("curve only cools", all(b[1] <= a[1] + 1e-9 for a, b in zip(curve, curve[1:]))),
        ("curve settles towards the surroundings",
         bool(curve) and ENV.temperature_k < curve[-1][1] < peak),
        ("per-second table has one row per second",
         _curve_rows(curve_table(result, 1.0, 5.0)) == 6),
    ]


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    start = time.perf_counter()
    failures: list[str] = []
    sources: dict[str, int] = {}

    print(f"{'result':>6}  input")
    print("-" * 78)
    for text, expected, source in CASES:
        result = analyse(text, ENV, 1.0)
        label = result.source or result.error or "?"
        sources[label] = sources.get(label, 0) + 1
        if result.error:
            failures.append(f"{text}: unexpected error: {result.error}")
            print(f"{'FAIL':>6}  {text}  ({result.error})")
            continue
        ok = result.equation == expected
        if source and source.lower() not in result.source.lower():
            ok = False
        if abs(result.mass_in - result.mass_out) > 1e-6 * max(1.0, result.mass_in):
            failures.append(
                f"{text}: mass {result.mass_in} in vs {result.mass_out} out"
            )
            ok = False
        for row in result.product_rows:
            if row["moles"] < -1e-12 or row["mass"] < -1e-12:
                failures.append(f"{text}: negative amount for {row['species'].display}")
                ok = False
        if result.final_ph is not None and result.final_ph.ph is not None:
            if not -5.0 <= result.final_ph.ph <= 20.0:
                failures.append(f"{text}: absurd pH {result.final_ph.ph}")
                ok = False
        print(f"{'ok' if ok else 'FAIL':>6}  {text}")
        if not ok:
            failures.append(
                f"{text}\n        expected: {expected}\n        got:      {result.equation}"
            )

    print("-" * 78)
    for text in ERRORS:
        result = analyse(text, ENV, 1.0)
        ok = result.error is not None
        print(f"{'ok' if ok else 'FAIL':>6}  {text or '<empty>'} -> {result.error or result.equation}")
        if not ok:
            failures.append(f"{text}: expected a clean error, got {result.equation}")

    print("-" * 78)
    checks = _redox_checks() + _temperature_checks() + _model_checks()
    for name, ok in checks:
        print(f"{'ok' if ok else 'FAIL':>6}  {name}")
        if not ok:
            failures.append(name)

    total = len(CASES) + len(ERRORS) + len(checks)
    print()
    print(f"cases: {total}   passed: {total - len(failures)}   failed: {len(failures)}")
    print("prediction sources used:")
    for name, count in sorted(sources.items(), key=lambda kv: -kv[1]):
        print(f"  {count:3d}  {name}")
    for item in failures:
        print("  !!", item)
    print(f"elapsed: {(time.perf_counter() - start) * 1000:.1f} ms")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
