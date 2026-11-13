"""MasChem: a terminal chemistry solver."""

from .solver import Analysis, analyse
from .surroundings import PRESETS, Surroundings, get

__version__ = "1.0.0"
__all__ = ["Analysis", "analyse", "PRESETS", "Surroundings", "get", "__version__"]
