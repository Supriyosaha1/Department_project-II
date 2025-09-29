"""
Rascas Analysis Package for Ramses and Gadget Output Processing

This package provides tools for analyzing astrophysical simulation output
from Ramses and Gadget codes using Rascas analysis techniques.
"""

__version__ = "0.1.0"
__author__ = "Department Project Team"

from .ramses import RamsesDataProcessor
from .gadget import GadgetDataProcessor
from .utils import ConfigManager, DataValidator

__all__ = [
    "RamsesDataProcessor",
    "GadgetDataProcessor", 
    "ConfigManager",
    "DataValidator"
]