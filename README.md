# Department Project II - Rascas Analysis for Ramses and Gadget

This repository provides a comprehensive Python framework for processing and analyzing astrophysical simulation output from **Ramses** and **Gadget** codes using **Rascas** analysis tools.

## Overview

The project enables seamless workflows for:
- Processing Ramses AMR simulation output
- Processing Gadget N-body/SPH simulation output  
- Preparing data for Rascas analysis
- Running automated Rascas analysis pipelines
- Visualizing and analyzing results

## Features

- **Modular Design**: Separate processors for Ramses and Gadget data
- **Flexible Configuration**: JSON and INI-based configuration management
- **Automated Workflows**: Command-line scripts for batch processing
- **Data Validation**: Built-in validation for simulation output
- **Format Conversion**: Convert between different data formats (e.g., Gadget to HDF5)
- **Rascas Integration**: Direct interface to Rascas analysis tools
- **Visualization**: Built-in plotting utilities for analysis results

## Installation

1. Clone the repository:
```bash
git clone https://github.com/Supriyosaha1/Department_project-II.git
cd Department_project-II
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Install the package:
```bash
pip install -e .
```

## Quick Start

### Processing Ramses Output

```bash
# Process all snapshots in a Ramses output directory
python scripts/process_ramses.py /path/to/ramses/output --output ./rascas_configs

# Process specific snapshots and run Rascas analysis
python scripts/process_ramses.py /path/to/ramses/output --snapshots 00001 00005 --run-rascas
```

### Processing Gadget Output

```bash
# Process Gadget snapshots with HDF5 conversion
python scripts/process_gadget.py /path/to/gadget/output --convert-hdf5 --output ./rascas_configs

# Process and run Rascas analysis
python scripts/process_gadget.py /path/to/gadget/output --run-rascas
```

### Python API Usage

```python
from rascas_analysis import RamsesDataProcessor, GadgetDataProcessor
from rascas_analysis.utils import ConfigManager, RascasRunner

# Load configuration
config = ConfigManager("config/default.json")

# Process Ramses data
ramses_processor = RamsesDataProcessor("/path/to/ramses/output")
snapshots = ramses_processor.discover_snapshots()
snapshot_data = ramses_processor.load_snapshot(snapshots[0])

# Prepare for Rascas analysis
config_file = ramses_processor.prepare_rascas_input(snapshot_data, "./output")

# Run Rascas analysis
runner = RascasRunner(config_manager=config)
result = runner.run_analysis(config_file)
```

## Directory Structure

```
Department_project-II/
├── rascas_analysis/           # Main package
│   ├── ramses/               # Ramses data processing
│   ├── gadget/               # Gadget data processing
│   └── utils/                # Utility modules
├── scripts/                  # Command-line scripts
├── config/                   # Configuration files
├── examples/                 # Example workflows
└── data/                     # Example data (create as needed)
```

## Configuration

The system uses flexible configuration management supporting both JSON and INI formats. Key configuration sections:

- **general**: Global settings (debug, parallel processing, etc.)
- **ramses**: Ramses-specific options (unit systems, loading preferences)
- **gadget**: Gadget-specific options (particle types, formats)
- **rascas**: Rascas executable and analysis settings
- **analysis**: Analysis parameters (binning, computations)
- **plotting**: Visualization settings

See `config/default.json` for a complete example configuration.

## Requirements

### Python Dependencies
- numpy >= 1.21.0
- scipy >= 1.7.0
- matplotlib >= 3.5.0
- astropy >= 5.0
- h5py >= 3.0.0
- pandas >= 1.3.0
- pynbody >= 1.0.0
- yt >= 4.0.0

### External Tools
- **Rascas**: The Rascas analysis executable must be installed and accessible in your PATH, or specify the path in configuration.

## Examples

Run the example workflow to test the installation:

```bash
python examples/example_workflow.py
```

This will demonstrate:
- Loading and validating simulation data
- Creating Rascas configurations
- Running analysis workflows
- Basic plotting and visualization

## Data Formats Supported

### Ramses
- AMR grid data (amr_*.out files)
- Particle data (part_*.out files)  
- Info files (info_*.txt)
- Multiple output directories (output_NNNNN)

### Gadget
- Binary snapshot files (snapshot_NNN)
- Header information parsing
- Multi-file snapshots
- HDF5 conversion support

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is developed for academic and research purposes.

## Support

For questions and support, please open an issue in the GitHub repository.

---

**Note**: This is a research tool designed for astrophysical simulation analysis. Ensure you have proper Rascas installation and simulation data before use.