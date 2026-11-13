# MasChem — usage reference

Everything you can type, and what it means. For the underlying model see
[how-it-works.md](how-it-works.md); for a five-minute introduction see the
[README](../README.md).

---

## 1. Running it

```
python -m maschem                          interactive shell
python -m maschem "2g Fe + 3g S"           one-shot answer
python -m maschem -e M1 "2g Fe + 3g S"     one-shot, on Mars
python -m maschem --selftest               built-in checks
python -m maschem M1                       print that environment
python -m maschem --help                   command list
```

Run from the folder that contains the `maschem` package. On Windows `maschem.bat` does the
`python -m maschem` part for you, so `maschem.bat "2g Fe + 3g S"` is equivalent.

| Flag | Meaning |
|---|---|
| `-e`, `--env`, `--surroundings <code or spec>` | start in that environment |
| `-v`, `--volume <litres>` | solution volume used for pH and thermal mass (default 1.0) |
| `--cooling`, `--loss <W/K>` | heat-loss coefficient for the cooling model (default 0.5) |
| `--no-colour` | no ANSI colour; `NO_COLOR=1` in the environment does the same |
| `--selftest` | run the regression checks and exit |
| `-h`, `--help` | print the command list |

Exit status: `0` success, `1` the reaction could not be analysed, `2` a bad flag or value.
Colour is only used when stdout is a terminal, so piping to a file or a pager is clean.

The report adapts to the terminal width (46–100 columns): long equations and atmosphere
lists wrap with aligned continuations, and below about 66 columns the atom column of the
product table is dropped.

## 2. What you can type

### 2.1 Two kinds of input

**A list of chemicals** — MasChem works out the products:

```
2g Fe + 3g S
10g Mg + 10g O2
KMnO4 + HCl
```

**A full equation** — you supply the products and MasChem checks and balances them:

```
CH4 + 2O2 -> CO2 + 2H2O
25g CaCO3 -> CaO + CO2
Fe + CuSO4 -> FeSO4 + Cu
```

Arrows may be written `->`, `-->`, `=>`, `→` or a spaced ` = `. Anything on the left is a
reactant, anything on the right a product.

### 2.2 Amounts and units

A number with a unit in front of a formula is an amount:

```
2g Fe          2 grams
25mg Zn        25 milligrams
1kg Al         1 kilogram
0.5 mol HCl    half a mole
250 mmol NaOH  250 millimoles
```

The words `grams`, `milligrams`, `kilograms`, `moles` and `millimoles` work as well as the
symbols, and `of` is optional: `2g of iron`, `2 grams of iron`, `2 g Fe`.

- **A bare number** in front of a formula is a stoichiometric coefficient, not an amount:
  `CH4 + 2O2 -> CO2 + 2H2O`.
- **No amounts at all** means the answer is per mole of the balanced reaction. The report
  says so in its notes.
- **Mixed input** is fine: given amounts decide the limiting reagent; species without
  amounts are consumed as needed.

Separators: `+`, and the words `mixed with`, `together with`, `reacts with`, `plus`, `and`.

### 2.3 Names, formulas, charges, phases

Common names are mapped to formulas before anything else: `water`, `table salt`, `iron`,
`copper`, `sulfur`/`sulphur`, `methane`, `ethanol`, `glucose`, `ammonia`, `carbon dioxide`
and about forty more (full list in `maschem/formula.py`, `NAME_ALIASES`).

Formulas support:

| Written | Parsed as |
|---|---|
| `Al2(SO4)3` | nested brackets with multipliers |
| `CuSO4·5H2O` | hydrates (the `·` or `.` separator with a leading count) |
| `NH4+`, `Fe3+`, `SO4^2-` | charges, with or without the caret |
| `MnCl2(aq)`, `FeS(s)` | phase labels `(s) (l) (g) (aq)` |

If you leave the phase out, MasChem infers it from the environment (solubility rules,
melting and boiling points at the current pressure) — that is why `Br2` is a liquid at
room temperature and `H2S` is a gas.

### 2.4 Environments on the line

Append `@` and a code to the end of any line to run *that line* somewhere else:

```
10g Mg + 10g O2 @S1
2g Fe + 3g S @M1
CH4 + 2O2 @T=350K P=0.5atm O2=100%
```

The session environment (set with the `surroundings` command) is used when no `@` suffix
is present.

## 3. Environments

### 3.1 Built-in codes

`surroundings list` prints this table live.

| Code | Place | T | Pressure | Atmosphere / water |
|---|---|---|---|---|
| `E1` | Earth surface | 288.15 K | 101.325 kPa | air, liquid water |
| `E2` | Earth orbit | 255 K | vacuum | — |
| `S1` | deep space | 2.73 K | vacuum | — |
| `ME1` | Mercury, day side | 440 K | vacuum | — |
| `V1` | Venus surface | 737 K | 9200 kPa | CO2 96.5 %, N2 3.5 % |
| `MO1` | Moon, day/night average | 250 K | vacuum | — |
| `M1` | Mars surface | 210 K | 610 Pa | CO2 95.1 %, N2 2.8 %, Ar 2 % |
| `J1` | Jupiter, 1 bar level | 165 K | 100 kPa | H2/He |
| `SA1` | Saturn, 1 bar level | 134 K | 100 kPa | H2/He |
| `T1` | Titan surface | 94 K | 146.7 kPa | N2/CH4 |
| `U1` | Uranus, 1 bar level | 76 K | 100 kPa | H2/He/CH4 |
| `N1` | Neptune, 1 bar level | 72 K | 100 kPa | H2/He/CH4 |
| `P1` | Pluto surface | 44 K | ~1 Pa | N2/CH4 |

Water is only liquid where the preset says so; anywhere else pH is undefined, ionic
products stay solid, and the report tells you.

### 3.2 Custom environments

```
surroundings T=350K P=0.5atm O2=100% water=off
surroundings T=90C P=2bar name=autoclave
surroundings P=0.001atm water=off
```

| Key | Meaning |
|---|---|
| `T=`, `temp=`, `temperature=` | temperature; plain number or `K` suffix is kelvin, `C` suffix is Celsius |
| `P=`, `pressure=` | pressure in `Pa`, `kPa`, `MPa`, `atm`, `bar`, `mbar`, `torr`, `mmHg` or `psi` (plain number = Pa) |
| `O2=`, `N2=`, `CO2=`, `Ar=`, `He=`, `H2=`, `CH4=`, `NH3=`, `SO2=` | gas fractions, as `%` or a fraction of 1 |
| `water=`, `h2o=` | `on`/`yes`/`true`/`1`/`liquid`/`available` or anything else for off |
| `name=`, `label=` | a label for the banner (underscores become spaces) |

Gases you name define the mixture; the leftover background gases share whatever fraction
is left, so `O2=100%` really is pure oxygen. Fractions adding up to more than 100 % are
rejected.

The `surroundings` command with no argument shows the current environment in full,
including where water would boil at that pressure and what neutral pH is.

## 4. Shell commands

| Command | Detail |
|---|---|
| `help`, `?`, `h` | command list and syntax |
| `examples` | print ready-made inputs |
| `surroundings` | show the current environment |
| `surroundings list` | table of every built-in environment |
| `surroundings <code>` | switch preset (`E1`, `M1`, …) |
| `surroundings <key=value …>` | switch to a custom environment, starting from the current one |
| `volume <litres>` | solution volume; sets pH concentration **and** thermal mass |
| `cooling` | show the current heat-loss coefficient |
| `cooling <W/K>` | set it; must be positive (0.05 dewar, 0.5 open beaker, 5 fan) |
| `dhf <species>` | read back the formation enthalpy in use |
| `dhf <species> <kJ/mol>` | add or override one for this session |
| `curve` | temperature and pH over the automatic window (4τ) |
| `curve <step> <end>` | the same at any resolution, in seconds |
| `clear` / `cls` | clear the screen |
| `quit` / `exit` / `q` | leave (Ctrl-C and Ctrl-D also work) |

`dhf` accepts the same formula syntax as input, phases included:

```
maschem> dhf MnCl2(aq)
MnCl2(aq): -555.1 kJ/mol
maschem> dhf MnSO4(aq) -1100
MnSO4(aq): -1100.0 kJ/mol  (this session only)
```

That value is then used for the rest of the session, which is how you get an energy and a
temperature out of a compound whose data was missing.

## 5. Reading the report

| Section | What it holds |
|---|---|
| surroundings | temperature, pressure, atmosphere, water phase at that pressure |
| reaction | your input, where the products came from, the balanced equation, limiting reagent |
| energy & temperature | ΔH per mole of reaction (negative = exothermic), total heat for the amounts given, the start temperature, the adiabatic **peak** — or **dip** if heat is absorbed — and the settled temperature |
| temperature & pH vs time | the cooling model, the window, one row per sample with t, T in K and °C, and pH; then ASCII graphs |
| pH | pH before and after mixing, the neutral value at that temperature, what drove it |
| colour | solution colour (from the ions present) and product colours |
| redox | per species: element, oxidation state before → after, electrons per atom, atoms involved, and which species was oxidised/reduced; electron total; standard potentials where known |
| products | species, phase, moles, mass |
| atoms | atoms and moles and mass of each element in the products (conservation check) |
| notes / warnings | assumptions, missing data, hazards |

Two conventions worth knowing: ΔH is **per mole of the reaction as balanced** (so
`2 KMnO4 + 16 HCl → …` is per two moles of permanganate), and the timing is "after the
reaction has released its heat", because the model has no kinetics.

Graphs use `·` for a sampled point, `│` to join consecutive samples, and the labelled rows
are the top, middle and bottom of the value range. The time scale is printed under the
axis.

## 6. Tuning the model

| Knob | Effect |
|---|---|
| `volume 0.025` | a test tube: pH concentrations rise, and the same energy now moves the temperature far more (and cools fast) |
| `cooling 0.05` | a dewar: τ grows 10×, the curve window grows with it |
| `cooling 5` | fan-cooled or thin-walled: τ shrinks 10× |
| `dhf <sp> <value>` | fills a thermochemistry gap for the session |
| `data.py` | the permanent place to add enthalpies, heat capacities, melting/boiling points, Ka values, colours and standard potentials |

Because the automatic window is four time constants, changing the cooling coefficient or
the volume changes *how much time* the table shows, not just the numbers — a 1 L beaker is
tracked for hours, a test tube for a minute.

## 7. Errors and what to do about them

| Message | Meaning |
|---|---|
| `unknown element 'Xx'` | the formula has an element that does not exist |
| `could not balance this equation` | no positive integer solution exists — check the formulas and the products |
| `a reactant would be left out` | the only balance zeroes a reactant; that is not the reaction you meant |
| `cannot work out the products` / `no reaction` | the rules do not justify a reaction; supply the products yourself |
| `no formation enthalpy for: …` | thermochemistry is missing; use `dhf`, or accept no temperature |
| `unknown surroundings 'X'` | not a preset code; use `surroundings list` or `T=…` |
| `unknown setting 'foo'` | custom environment keys are limited to those in §3.2 |

Deliberate refusals (for example `Al2O3 + NaOH`, amphoteric; `ZnO + H2`, no reduction)
are a feature: MasChem would rather say "I cannot justify this" than print a forecast.

## 8. Recipes

```
surroundings list                                  see everything available
surroundings M1                                    work on Mars from now on
10g Mg + 10g O2                                    hot solid reaction
2g Fe + 3g S @M1                                   a one-off on Mars
volume 0.05                                        a 50 mL flask
cooling 2                                          open bench, some draught
KMnO4 + HCl                                        dip, then recovery, with pH
curve 5 3600                                       that curve every 5 s for an hour
dhf MnSO4(aq) -1100                                fill a data gap
2 mol Mn + 2 mol H2SO4 -> 2 mol MnSO4 + 2 mol H2    ...now it has a temperature
surroundings T=500K water=off                      dry, hot, no solvent
```
