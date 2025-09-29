#!/usr/bin/env python3
"""
Example workflow for Rascas analysis with Ramses and Gadget data

This script demonstrates how to use the rascas_analysis package to process
simulation output and run Rascas analysis.
"""

import sys
from pathlib import Path
import logging

# Add the parent directory to the path so we can import our package
sys.path.insert(0, str(Path(__file__).parent.parent))

from rascas_analysis import RamsesDataProcessor, GadgetDataProcessor
from rascas_analysis.utils import ConfigManager, RascasRunner, DataValidator


def main():
    """Main function to run examples."""
    print("Rascas Analysis Example Workflow")
    print("================================")
    print("This is a demonstration workflow for processing Ramses and Gadget data.")
    print("To use with real data, update the data paths in the script.")


if __name__ == '__main__':
    main()
