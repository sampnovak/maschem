# How MasChem works

This is the model behind the numbers, in the order the program executes it. Nothing here
is looked up as a finished answer: the data tables hold *physical constants* (enthalpies,
heat capacities, melting points, Ka values, potentials), and everything else — products,
coefficients, temperatures, pH, electrons — is computed from them.

```
text ─▶ parse ─▶ predict products ─▶ balance ─▶ infer phases
                                            │
              ┌─────────────────────────────┴───────────────────────────┐
              ▼                                                         ▼
      thermochemistry                                            aqueous chemistry
  ΔH, adiabatic peak/dip                                pH before/after, hydrolysis
              │                                                         │
              ▼                                                         ▼
      Newton cooling  ───────▶  time series: T(t), pH(t)  ◀─────  pH at each T
              │
              ▼
      redox · accounting · colour · report
```

## 1. Parsing

`formula.py` turns what you typed into a `Species`: element counts, ionic charge and an
optional phase. Brackets are multiplied out, hydrates (`CuSO4·5H2O`) are expanded into
their water, charges (`NH4+`, `SO4^2-`) are separated from the formula, common names are
replaced by formulas before parsing, and a phase label overrides everything.

`solver.py` splits the line on `+` (and the words *and*, *plus*, *mixed with*, …) and then
on an arrow. A leading number with a mass or mole unit becomes an amount; a bare number is
a stoichiometric coefficient.

## 2. Extent of reaction and the limiting reagent

For amounts `n_i` and coefficients `ν_i`, the extent is

```
ξ = min over reactants of (n_i / ν_i)
```

The reactant that minimises it is reported as limiting. Products get `n = ν·ξ`, and any
reactant with `n_i − ν_i·ξ > 0` is reported as left over. With no amounts at all, `ξ = 1`
and everything is per mole of the balanced reaction.

## 3. Predicting products

The database is consulted first: if the exact reactant set has a known reaction, its
products are used. If not, rules are tried in a fixed order, cheapest first:

1. combustion (hydrocarbon + O2 → CO2 + H2O, assumed complete)
2. metal + non-metal → binary ionic compound
3. element + O2 → oxide
4. acid + base → salt + water
5. acid + ammonia/amine → ammonium salt
6. acidified oxidiser (KMnO4, K2Cr2O7, MnO2, PbO2 with HCl/HBr/HI) → halide oxidised
7. metal oxide + acid → salt + water
8. acidic oxide + base → salt + water
9. acid + sulfite → salt + water + SO2
10. acid + sulfide → salt + H2S
11. metal oxide + H2 → metal + water (only where reduction is possible)
12. acid + carbonate/bicarbonate → salt + water + CO2
13. metal + acid → salt + H2 (only above hydrogen in the activity series)
14. oxide + water → hydroxide or acid
15. double displacement (salt + salt, salt + hydroxide) → precipitate if insoluble
16. single displacement (metal + salt) → the more reactive metal wins
17. carbonate decomposition → oxide + CO2
18. hydroxide decomposition → oxide + water

Each rule is a *generalisation*, not a list of examples. Rule 12 asks "is this a salt whose
anion is CO3 or HCO3?", so `Li2CO3`, `(NH4)2CO3` and a carbonate invented on the spot are
all handled. Rule 6 knows the half-reactions, so the oxidation state of the halogen comes
out right without a table of finished equations.

When no rule applies the answer is an error, never a guess: `Al2O3 + NaOH` (amphoteric)
and `ZnO + H2` (no reduction) are refused on purpose.

## 4. Balancing

The equation is a null space problem. Build a matrix with one row per element plus one row
for charge, and one column per species (reactants negative, products positive). Balancing
means finding a positive integer vector `ν` with `A·ν = 0`.

- The arithmetic is exact: every entry is a `Fraction`, reduced row echelon form is
  computed with rationals, so no floating-point drift can fake a solution.
- The null-space basis is turned into the smallest whole numbers (clear denominators, then
  divide by the GCD).
- Solutions that zero out a species you wrote are rejected as degenerate — that is the
  `a reactant would be left out` error.
- Charge is a row, so ionic equations balance in charge as well as in atoms.

## 5. Phases

With no phase label, `infer_state` decides: water follows the triple point and the boiling
point at the current pressure; dissolved ions and soluble salts are `aq` when there is
liquid water; strong/weak acids and bases are solutions; volatile ones (H2S, NH3) become
gases when they are hot; everything else uses tabulated melting and boiling points.

Boiling points off 1 atm come from Clausius–Clapeyron,

```
1/T_bp(P) = 1/T_bp(1 atm) − R·ln(P/1 atm)/ΔH_vap
```

with ΔH_vap from Trouton's rule (≈ 88·T_bp) when it is not tabulated. That is why `E2`,
`M1` and `S1` give different phases from `E1`, and why water boils at 373.1 K on Earth but
not at all on Mars.

## 6. Thermochemistry

### Enthalpy of reaction

```
ΔH = Σ_products ν·ΔHf° − Σ_reactants ν·ΔHf°
```

ΔHf° comes from `data.py`, keyed by formula *and* phase (`MnCl2(s)`, `MnCl2(aq)`), so
dissolving a salt changes the energy. Elements in their standard state are 0 by
definition. If any species is missing, the reaction gets no enthalpy — and a warning that
names the missing species — instead of a wrong number.

### Adiabatic temperature, solid/gas reactions

Heat released raises the products' temperature. Heat capacity is tabulated where known and
estimated from atom count where not (25 + 4·N J/mol·K for gases, ~25·N for solids). Above a
product's boiling point the model switches to an approximate high-temperature scaling,

```
Cp(T) ≈ Cp° · (T/298.15)^0.25,
```

and the temperature is found by iterating

```
T_peak = T_0 − ΔH·1000 / Cp(T_peak)
```

with damping until the change is below 0.5 K. If the fixed point falls below absolute zero
the reaction cannot be self-sustaining and the report says so. The result is an **upper
bound**: no dissociation, no latent heat, no heat loss.

### Aqueous reactions

When water is the medium, the solution — not the stoichiometric product mixture — is what
absorbs the heat:

```
ΔT = −ΔH_total·1000 / (4180 · V)          T_peak = T_env + ΔT
```

where 4180 J/(L·K) is 1 L of water at 4.18 J/(g·K). Because the sign is kept, an
endothermic reaction produces a **dip**: permanganate in hydrochloric acid goes from 288 K
down to 261 K and then warms back up. That is physical (the reaction is entropy-driven),
and the report labels it a dip rather than pretending a peak exists.

## 7. Cooling

The reaction is treated as instantaneous, so the mixture starts at its peak (or dip) and
relaxes towards the surroundings:

```
C_total = 4180·V              (aqueous)   or   Σ ν·Cp(products)·ξ   (dry)
τ       = C_total / h
T(t)    = T_env + (T_peak − T_env)·e^(−t/τ)
```

`h` is the heat-loss coefficient in W/K: 0.5 by default, `cooling 0.05` for a dewar,
`cooling 5` for something fan-cooled. τ is the time to close 63 % of the gap; the
automatic table runs for **4τ** (98 % recovered) and picks the coarsest "nice" step
(1, 2, 5, 10, 15, 30 s, 1, 2, 5, 10, 15, 30 min …) that still fits the window in about
24 rows. That is why a litre of water is tracked for hours and a test tube for seconds,
and why the table never shows twenty identical numbers again: the window scales with the
vessel, not with a fixed 60 seconds.

Some consequences worth remembering:

- Doubling `volume` doubles `C_total`, so τ doubles and the same energy gives half the ΔT.
- A bigger `cooling` value shortens τ but does not change the peak.
- In vacuum the model keeps using `h` and warns that radiation, not convection, is the
  real mechanism.

## 8. pH

pH comes from a charge balance on the final mixture. With `[H+] = h` and `[OH−] = Kw/h`:

```
h + [BH+]  =  Kw/h + [A−] + [strong acid]
```

- strong acid: `[A−] = C`
- weak acid (Ka): `[A−] = C·Ka/(Ka + h)`
- weak base or hydrolysing anion (Kb): `[BH+] = C·h/(h + Kw/Kb)`
- a dissolved salt contributes the hydrolysis of its own ions (NH4+, CO3²⁻, CH3COO⁻, …)

The equation is solved for log10 h by bisection over 10⁻¹⁶ … 10⁻² M — no approximations
about which term dominates. Concentrations come from moles over the solution volume, so
`volume` changes pH as well as thermal mass.

Kw is temperature dependent (van 't Hoff, ΔH_ionisation = 55.8 kJ/mol):

```
Kw(T) = Kw(298) · exp( −(ΔH/R)·(1/T − 1/298.15) ),   neutral pH = pKw(T)/2
```

which is what makes the pH curve move: the mixture is fully reacted at t = 0, so pH(t)
changes only because Kw (and the conjugate-base constants derived from it) change with the
temperature of the cooling solution. The report says this explicitly under the graph.

## 9. Redox

Oxidation states are assigned per element inside each species: single-element species from
their charge, a small override table for the awkward cases (peroxides O −1, superoxide
O −½, mixed oxides like Fe3O4 → Fe 8/3, N2H4, NH4NO3, glucose), otherwise fluorine −1,
oxygen −2 (except H2O2), alkali metals +1, alkaline earths +2, hydrogen +1 (or −1 in
hydrides), and the remainder follows from the species charge.

Then the same element is matched **species by species** between reactants and products:

- atoms whose state does not change are cancelled first (chloride in `2 KMnO4 + 16 HCl`:
  only ten of the sixteen chlorines are oxidised),
- the rest are paired greedily, giving `from_species → to_species` for each change,
- electrons per atom are `|state_after − state_before|`, and the total is that times the
  atoms moved.

Standard potentials from a 32-couple table (`E0`) give `E°cell = E°cathode − E°anode` for
the noted couples, which is used as a feasibility cross-check on the rule-based
predictions.

## 10. Accounting, colour, layout

Mass is conserved by construction (the balance does it); the report verifies it as
`mass in = mass out`. Atom counts are `moles × atoms-per-formula × 6.02214076×10²³`.
Colour comes from the ions present (Cu²⁺ blue, MnO4⁻ deep purple, Fe³⁺ yellow, …) with a
fallback to compound and gas colours, so a salt solution is coloured by what is actually
dissolved.

Rendering measures the terminal once per report and lays out to that width, wrapping long
values with aligned continuations. Graphs are plain ASCII: a column per time slice,
interpolated between samples, with `·` for sampled points and `│` joining them.

## 11. Where the approximations are

| Model | Assumption | Consequence |
|---|---|---|
| Kinetics | instantaneous reaction | no rates, no induction period, no "half-way" states |
| Adiabatic peak | no loss, no dissociation, no latent heat | an upper bound; the warning appears above 3000 K |
| Heat capacity | tabulated, else atom-count estimate | estimated Cp is listed in the notes |
| Cooling | single lumped coefficient `h` | no geometry, no conduction gradients, no radiation detail |
| pH | concentrations, no activity coefficients | good to ~0.1 pH in dilute solution |
| Solubility | qualitative rules | `sparingly soluble` cases are simplified |
| Redox | standard potentials only | no Nernst correction for concentration |

## 12. Extending it

- **A missing enthalpy** — `dhf <species> <kJ/mol>` for a session, or add it to `DHF` in
  `data.py` (key with phase: `MnCl2(aq)`). Add `CP`, `MP`, `BP` the same way.
- **A missing rule** — add a block to `predict()` in `reactions.py`; return a
  `Prediction(products, "your rule", notes)`. Keep the database lookup ahead of it.
- **A new body** — add an entry to `_BODY_TABLE` in `surroundings.py`; `PRESETS` and
  `surroundings list` pick it up automatically.
- **A check for it** — add a case to `selftest.py` so the behaviour stays.
