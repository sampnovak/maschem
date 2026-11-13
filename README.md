# MasChem

A chemistry solver that lives in your terminal. Type in the chemicals and you get the
balanced equation, the heat involved, how hot (or cold) the mixture actually gets, how
that temperature decays second by second, the pH over the same timeline, the electrons
that moved, what is left over and what colour it all is — for surroundings anywhere from
Earth's surface to Pluto.

Pure Python, standard library only. No install, no accounts, no network.

```
maschem> 2g Fe + 3g S
maschem> 10g Mg + 10g O2 @S1
maschem> KMnO4 + HCl
```

---

## What it does

- **Reads what you would write.** `2g Fe + 3g S`, `0.5 mol HCl + 0.5 mol NaOH`,
  `2g of iron mixed with 2g of sulfur`, or a full equation with an arrow.
- **Predicts the products** when you leave them out — a reaction database first, then
  about eighteen inorganic rules (combustion, neutralisation, precipitation,
  displacement, acidified oxidisers, carbonate/sulfide/sulfite reactions,
  decompositions, …). Nothing is looked up by name: a carbonate is recognised as
  *a salt whose anion is CO3*, so `Li2CO3 + HCl` works even though nobody typed it in.
- **Balances exactly** with integer linear algebra (rational arithmetic, then the
  smallest whole numbers). It refuses nonsense instead of inventing a solution.
- **Works out the energy** from formation enthalpies and heat capacities, including the
  adiabatic peak — or the *dip*, when the reaction absorbs heat.
- **Shows temperature over time** with a Newton-cooling model whose window scales to the
  vessel: a 1 L beaker of permanganate is tracked for 9 hours, a gram of thermite for
  seconds. A `curve` command re-prints the same table at any resolution you like.
- **Shows pH over the same timeline**, recomputed at every temperature, plus solution and
  solid colours.
- **Explains the redox** per species: every oxidation state, per-element electrons, the
  electron total, and standard potentials where they are known.
- **Accounts for everything**: limiting reagent, leftovers, mass in = mass out, atom
  counts, and the phase of each product — including at the peak temperature.
- **Runs anywhere**: Earth's surface, Earth orbit, deep space, Mercury, Venus, the Moon,
  Mars, Jupiter, Saturn, Titan, Uranus, Neptune, Pluto, or an atmosphere you invent.

## Requirements

- Python 3.10 or newer (developed on 3.14). No third-party packages.
- Run it from the folder that *contains* the `maschem` package:

```
cd G:\glassPickle\MasChem\MasChem
python -m maschem
```

On Windows you can also just run `maschem.bat` in that folder. Any of the arguments below
can be appended to either form.

## Quick start

| You type | You get |
|---|---|
| `python -m maschem` | the interactive shell (`maschem>`) |
| `python -m maschem "2g Fe + 3g S"` | one answer, then exit |
| `python -m maschem -e M1 "2g Fe + 3g S"` | the same reaction on Mars |
| `python -m maschem --selftest` | the built-in solver checks |
| `python -m maschem --help` | the command list |

Once inside the shell, `help` prints the same command list, `examples` prints things worth
typing, and `quit` (or Ctrl-C) leaves.

## A real result

```
$ python -m maschem "2g Fe + 3g S"

MAS CHEM  ·  terminal reaction solver
surroundings E1 · Earth surface (sea level, room temperature)

─ SURROUNDINGS ─────────────────────────────────────────────
  temperature   288.15 K  (15.0 °C)
  pressure      101.325 kPa  (1.000 atm)
  atmosphere    N2 78.08% · O2 20.95% · Ar 0.93% · CO2 0.04%
  water         liquid (boils at 373.1 K here)

─ REACTION ─────────────────────────────────────────────────
  input         2g Fe + 3g S
  source        reaction database
  balanced      Fe(s) + S(s) → FeS(s)
  limiting      Fe(s)

─ ENERGY & TEMPERATURE ─────────────────────────────────────
  enthalpy      -100.0 kJ per mol of reaction  (exothermic)
  heat          3.58 kJ released
  t = 0         288.15 K  (surroundings)
  peak          2276 K  (adiabatic estimate, +1988 K)
  settled       288.15 K  (after heat exchange with surroundings)
  phases        FeS(s) → l at 2276 K

─ TEMPERATURE & pH VS TIME ─────────────────────────────────
  model         Newton cooling, τ = 3.6 s (1.8 J/K, 0.5 W/K loss to the room)
  window        14 s (4τ), sampled every 1 s
       t/s      T/K       °C
       0.0   2276.2   2003.1
       1.0   1794.4   1521.2
       2.0   1429.3   1156.2
       3.0   1152.7    879.6
       5.0    784.4    511.3
      10.0    412.0    138.9
      14.0    329.0     55.8
  temperature
  2276.2 ┤··│
         ┤  ··│
         ┤    ···│
         └────────────────────────────
          0 s                  14 s
```

(Shortened: the real report also has pH, colour, redox, products, atoms and notes.)

## Things to try

```
10g Mg + 10g O2                 thermite-hot solid reaction, per-gram amounts
KMnO4 + HCl                     acidified oxidiser, per-species redox, a dip in T
KMnO4 + HBr                     the same oxidiser, exothermic instead
2g of iron mixed with 2g of sulfur
0.5 mol HCl + 0.5 mol NaOH      neutralisation, pH before and after
25g CaCO3 -> CaO + CO2          you supply the products, MasChem checks them
Fe + CuSO4 -> FeSO4 + Cu        displacement
10g Mg + 10g O2 @S1             no air, no water, no convection
2g Fe + 3g S @M1                Mars pressure and temperature
curve 1 120                     last reaction, temperature and pH every second
surroundings list               every built-in environment
```

## Commands in the shell

| Command | What it does |
|---|---|
| `help` | command list and input syntax |
| `examples` | ready-made inputs |
| `surroundings` | show the current environment |
| `surroundings list` | table of every built-in environment |
| `surroundings M1` | switch environment (any code below) |
| `surroundings T=350K P=0.5atm O2=100% water=off` | build a custom one |
| `volume 2` | solution volume in litres — sets pH strength **and** thermal mass |
| `cooling 0.5` | heat loss to the room in W/K (0.05 dewar, 0.5 beaker, 5 fan) |
| `dhf MnCl2(aq) -555.1` | add a formation enthalpy for a compound MasChem lacks |
| `curve` | temperature and pH over the automatic window (4τ) |
| `curve 1 120` | …every 1 s for 120 s, with graphs |
| `clear` / `quit` | clear the screen / leave |

## Command-line flags

| Flag | Meaning |
|---|---|
| `-e`, `--env`, `--surroundings <code\|spec>` | start in that environment |
| `-v`, `--volume <litres>` | solution volume (default 1.0) |
| `--cooling`, `--loss <W/K>` | heat-loss coefficient (default 0.5) |
| `--no-colour` | plain text (also honours `NO_COLOR=1`) |
| `--selftest` | run the built-in checks |
| `-h`, `--help` | the command list |

Exit status is `0` on success, `1` if the reaction could not be analysed, `2` for a bad
flag. Any single environment code on its own (`python -m maschem M1`) prints that
environment instead of running a reaction.

## Environments

Append `@CODE` to a line to use one for that line only, or set it for the session with
`surroundings CODE`. `E1`, `E2` and `S1` are the ones you will use most; the rest of the
solar system comes free.

| Code | Place | Temperature | Pressure |
|---|---|---|---|
| `E1` | Earth surface | 288.15 K | 1 atm air, water |
| `E2` | Earth orbit | 255 K | vacuum |
| `S1` | deep space | 2.73 K | vacuum |
| `ME1` | Mercury (day) | 440 K | vacuum |
| `V1` | Venus surface | 737 K | 90.8 atm CO2 |
| `MO1` | Moon (average) | 250 K | vacuum |
| `M1` | Mars surface | 210 K | 0.006 atm CO2 |
| `J1` | Jupiter (1 bar) | 165 K | 1 bar H2/He |
| `SA1` | Saturn (1 bar) | 134 K | 1 bar |
| `T1` | Titan surface | 94 K | 1.45 atm N2/CH4 |
| `U1` | Uranus (1 bar) | 76 K | 1 bar |
| `N1` | Neptune (1 bar) | 72 K | 1 bar |
| `P1` | Pluto surface | 44 K | near-vacuum N2 |

Custom environments take `T=`, `P=`, `O2=`, `N2=`, `CO2=`, `Ar=`, `He=`, `H2=`, `CH4=`,
`NH3=`, `SO2=`, `water=on|off` and `name=`. Full syntax in
[docs/usage.md](docs/usage.md).

## Where the numbers come from

Enthalpies, heat capacities, melting/boiling points, acid constants, solubilities,
standard potentials and colours are data tables in `maschem/data.py`. Everything else is
computed: the balancing is linear algebra, the peak temperature is an energy balance, the
cooling curve is an exponential decay, pH is a charge balance solved numerically. If a
compound has no data, MasChem says so rather than guessing — and `dhf` lets you fill the
gap at runtime.

How each step works, with the equations, is in [docs/how-it-works.md](docs/how-it-works.md).

## What it does not do

Being honest about the model is part of using it well:

- **No kinetics.** The reaction is instantaneous; there is no rate, no induction period,
  no half-life. "Temperature at t = 5 s" means "after the heat was released", not "the
  reaction finished at 5 s".
- **The peak is an upper bound.** It is adiabatic: no heat escapes, nothing dissociates,
  nothing boils away. Real temperatures are lower.
- **No latent heat, no radiative detail.** Boiling and melting of the products are not
  cooled by enthalpy of vaporisation; in a vacuum the model still uses a conductive loss
  coefficient (with a warning that radiation is the real mechanism).
- **Limited thermochemistry.** About 220 species have formation enthalpies (and 110
  heat capacities). Missing ones are listed by name; `dhf` adds them for the session.
- **Rule-based predictions are inorganic-leaning.** Organic products beyond complete
  combustion (CO2 + H2O) are not predicted.
- **pH is a textbook model.** Strong/weak acids and bases, salt hydrolysis, one litre of
  water unless you say otherwise, no activity coefficients.
- **Deliberate refusals.** Reactions the rules cannot justify (amphoteric `Al2O3 + NaOH`,
  `ZnO + H2`) produce an error rather than a fabricated answer.

## Layout

This folder is the `maschem` package itself; `maschem.bat` sits one folder up, beside it,
and that parent folder is the one you run from.

```
README.md            this file
docs/
  usage.md           full input, command and environment reference
  how-it-works.md    the maths and the rules behind each number

__main__.py          entry point
cli.py               shell, commands, flags
solver.py            the pipeline: parse → predict → balance → energy → time series
formula.py           formula/name/charge/phase parsing
balance.py           exact equation balancing
reactions.py         product database and the inorganic rules
thermo.py            enthalpy, adiabatic temperature, Newton cooling
aqueous.py           pH and hydrolysis
redox.py             oxidation states, electrons, potentials
surroundings.py      environments (presets and custom)
data.py              the data tables
report.py            terminal layout, tables and graphs
selftest.py          built-in regression checks
```

## Tests

```
python -m maschem --selftest
```

Seventy-odd checks covering product predictions, error handling, balancing, per-species
redox, thermochemistry, the cooling curve, pH time series, layout and the environment
presets. It finishes in well under a second and exits non-zero if anything fails, so it
can be wired into a commit hook.
