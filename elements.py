"""Periodic-table data used by MasChem.

Atomic masses are standard atomic weights in u (bracket values are mass
numbers for elements with no stable isotope).  Electronegativity is the
Pauling scale (0.0 where no useful value exists).  ``states`` lists common
oxidation states, most common first.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Element:
    symbol: str
    name: str
    mass: float
    electronegativity: float
    states: tuple[int, ...]


_RAW: dict[str, tuple[str, float, float, tuple[int, ...]]] = {
    "H": ("Hydrogen", 1.008, 2.20, (1, -1)),
    "He": ("Helium", 4.0026, 0.0, ()),
    "Li": ("Lithium", 6.94, 0.98, (1,)),
    "Be": ("Beryllium", 9.0122, 1.57, (2,)),
    "B": ("Boron", 10.81, 2.04, (3,)),
    "C": ("Carbon", 12.011, 2.55, (4, -4, 2)),
    "N": ("Nitrogen", 14.007, 3.04, (-3, 3, 5, 4, 2)),
    "O": ("Oxygen", 15.999, 3.44, (-2,)),
    "F": ("Fluorine", 18.998, 3.98, (-1,)),
    "Ne": ("Neon", 20.180, 0.0, ()),
    "Na": ("Sodium", 22.990, 0.93, (1,)),
    "Mg": ("Magnesium", 24.305, 1.31, (2,)),
    "Al": ("Aluminium", 26.982, 1.61, (3,)),
    "Si": ("Silicon", 28.085, 1.90, (4, -4)),
    "P": ("Phosphorus", 30.974, 2.19, (5, 3, -3)),
    "S": ("Sulfur", 32.06, 2.58, (-2, 4, 6)),
    "Cl": ("Chlorine", 35.45, 3.16, (-1, 1, 3, 5, 7)),
    "Ar": ("Argon", 39.948, 0.0, ()),
    "K": ("Potassium", 39.098, 0.82, (1,)),
    "Ca": ("Calcium", 40.078, 1.00, (2,)),
    "Sc": ("Scandium", 44.956, 1.36, (3,)),
    "Ti": ("Titanium", 47.867, 1.54, (4, 3)),
    "V": ("Vanadium", 50.942, 1.63, (5, 4, 3)),
    "Cr": ("Chromium", 51.996, 1.66, (3, 6, 2)),
    "Mn": ("Manganese", 54.938, 1.55, (2, 4, 7)),
    "Fe": ("Iron", 55.845, 1.83, (2, 3)),
    "Co": ("Cobalt", 58.933, 1.88, (2, 3)),
    "Ni": ("Nickel", 58.693, 1.91, (2,)),
    "Cu": ("Copper", 63.546, 1.90, (2, 1)),
    "Zn": ("Zinc", 65.38, 1.65, (2,)),
    "Ga": ("Gallium", 69.723, 1.81, (3,)),
    "Ge": ("Germanium", 72.630, 2.01, (4,)),
    "As": ("Arsenic", 74.922, 2.18, (5, 3, -3)),
    "Se": ("Selenium", 78.971, 2.55, (-2, 4, 6)),
    "Br": ("Bromine", 79.904, 2.96, (-1, 1, 5)),
    "Kr": ("Krypton", 83.798, 0.0, ()),
    "Rb": ("Rubidium", 85.468, 0.82, (1,)),
    "Sr": ("Strontium", 87.62, 0.95, (2,)),
    "Y": ("Yttrium", 88.906, 1.22, (3,)),
    "Zr": ("Zirconium", 91.224, 1.33, (4,)),
    "Nb": ("Niobium", 92.906, 1.60, (5,)),
    "Mo": ("Molybdenum", 95.95, 2.16, (6, 4)),
    "Tc": ("Technetium", 98.0, 1.90, (7,)),
    "Ru": ("Ruthenium", 101.07, 2.20, (3, 4)),
    "Rh": ("Rhodium", 102.91, 2.28, (3,)),
    "Pd": ("Palladium", 106.42, 2.20, (2,)),
    "Ag": ("Silver", 107.87, 1.93, (1,)),
    "Cd": ("Cadmium", 112.41, 1.69, (2,)),
    "In": ("Indium", 114.82, 1.78, (3,)),
    "Sn": ("Tin", 118.71, 1.96, (4, 2)),
    "Sb": ("Antimony", 121.76, 2.05, (5, 3, -3)),
    "Te": ("Tellurium", 127.60, 2.10, (-2, 4, 6)),
    "I": ("Iodine", 126.90, 2.66, (-1, 1, 5, 7)),
    "Xe": ("Xenon", 131.29, 0.0, ()),
    "Cs": ("Caesium", 132.91, 0.79, (1,)),
    "Ba": ("Barium", 137.33, 0.89, (2,)),
    "La": ("Lanthanum", 138.91, 1.10, (3,)),
    "Ce": ("Cerium", 140.12, 1.12, (3, 4)),
    "Pr": ("Praseodymium", 140.91, 1.13, (3,)),
    "Nd": ("Neodymium", 144.24, 1.14, (3,)),
    "Pm": ("Promethium", 145.0, 1.13, (3,)),
    "Sm": ("Samarium", 150.36, 1.17, (3,)),
    "Eu": ("Europium", 151.96, 1.20, (3, 2)),
    "Gd": ("Gadolinium", 157.25, 1.20, (3,)),
    "Tb": ("Terbium", 158.93, 1.20, (3,)),
    "Dy": ("Dysprosium", 162.50, 1.22, (3,)),
    "Ho": ("Holmium", 164.93, 1.23, (3,)),
    "Er": ("Erbium", 167.26, 1.24, (3,)),
    "Tm": ("Thulium", 168.93, 1.25, (3,)),
    "Yb": ("Ytterbium", 173.05, 1.10, (3, 2)),
    "Lu": ("Lutetium", 174.97, 1.27, (3,)),
    "Hf": ("Hafnium", 178.49, 1.30, (4,)),
    "Ta": ("Tantalum", 180.95, 1.50, (5,)),
    "W": ("Tungsten", 183.84, 2.36, (6,)),
    "Re": ("Rhenium", 186.21, 1.90, (7,)),
    "Os": ("Osmium", 190.23, 2.20, (4,)),
    "Ir": ("Iridium", 192.22, 2.20, (3,)),
    "Pt": ("Platinum", 195.08, 2.28, (2, 4)),
    "Au": ("Gold", 196.97, 2.54, (3, 1)),
    "Hg": ("Mercury", 200.59, 2.00, (2, 1)),
    "Tl": ("Thallium", 204.38, 1.62, (1, 3)),
    "Pb": ("Lead", 207.2, 2.33, (2, 4)),
    "Bi": ("Bismuth", 208.98, 2.02, (3,)),
    "Po": ("Polonium", 209.0, 2.00, (4, 2)),
    "At": ("Astatine", 210.0, 2.20, (-1,)),
    "Rn": ("Radon", 222.0, 0.0, ()),
    "Fr": ("Francium", 223.0, 0.70, (1,)),
    "Ra": ("Radium", 226.0, 0.90, (2,)),
    "Ac": ("Actinium", 227.0, 1.10, (3,)),
    "Th": ("Thorium", 232.04, 1.30, (4,)),
    "Pa": ("Protactinium", 231.04, 1.50, (5,)),
    "U": ("Uranium", 238.03, 1.38, (6, 4)),
}

ELEMENTS: dict[str, Element] = {
    sym: Element(sym, name, mass, en, states)
    for sym, (name, mass, en, states) in _RAW.items()
}

METALS = frozenset(
    {
        "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
        "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
        "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce",
        "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
        "Bi", "Po", "Fr", "Ra", "Ac", "Th", "Pa", "U",
    }
)

NONMETALS = frozenset(
    {"H", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I", "Te", "At"}
)

NOBLE_GASES = frozenset({"He", "Ne", "Ar", "Kr", "Xe", "Rn"})

HALOGENS = frozenset({"F", "Cl", "Br", "I", "At"})

ALKALI = frozenset({"Li", "Na", "K", "Rb", "Cs", "Fr"})
ALKALINE = frozenset({"Be", "Mg", "Ca", "Sr", "Ba", "Ra"})

# Metals roughly in order of decreasing reactivity (activity series).
REACTIVITY_SERIES = [
    "Cs", "Fr", "Rb", "K", "Ba", "Sr", "Ca", "Na", "Li", "Mg", "Al", "Mn",
    "Zn", "Cr", "Fe", "Cd", "Co", "Ni", "Sn", "Pb", "H", "Cu", "Hg", "Ag",
    "Pd", "Pt", "Au",
]

DIATOMIC = {
    "H": "H2",
    "N": "N2",
    "O": "O2",
    "F": "F2",
    "Cl": "Cl2",
    "Br": "Br2",
    "I": "I2",
}


def element(symbol: str) -> Element:
    try:
        return ELEMENTS[symbol]
    except KeyError:
        raise KeyError(f"unknown element symbol {symbol!r}") from None


def is_metal(symbol: str) -> bool:
    return symbol in METALS
