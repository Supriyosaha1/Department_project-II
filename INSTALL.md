# Installation Guide

## Prerequisites

- Python 3.8 or higher
- pip package manager
- Git

## Quick Installation

1. **Clone the repository:**
```bash
git clone https://github.com/Supriyosaha1/Department_project-II.git
cd Department_project-II
```

2. **Install Python dependencies:**
```bash
pip install -r requirements.txt
```

3. **Install the package in development mode:**
```bash
pip install -e .
```

## Dependencies

### Core Scientific Libraries
- numpy >= 1.21.0
- scipy >= 1.7.0
- matplotlib >= 3.5.0
- astropy >= 5.0
- h5py >= 3.0.0

### Data Analysis
- pandas >= 1.3.0
- seaborn >= 0.11.0
- plotly >= 5.0.0

### Astrophysical Tools
- pynbody >= 1.0.0
- yt >= 4.0.0

### Utilities
- tqdm >= 4.62.0
- configparser >= 5.0.0

## External Requirements

### Rascas Installation

The Rascas analysis tool must be installed separately. Options include:

1. **Install from source** (recommended for custom builds)
2. **Use package manager** (if available for your system)
3. **Specify executable path** in configuration

#### Configuration

Update the `config/default.json` file to specify the Rascas executable path:

```json
{
  "rascas": {
    "executable": "/path/to/rascas",
    "config_template": "default.cfg",
    "output_format": "hdf5"
  }
}
```

## Verification

Test the installation:

```bash
# Test basic imports
python -c "from rascas_analysis import RamsesDataProcessor, GadgetDataProcessor; print('✓ Installation successful')"

# Run example workflow
python examples/example_workflow.py
```

## Troubleshooting

### Common Issues

1. **Import errors:** Ensure all dependencies are installed
2. **Rascas not found:** Update the executable path in configuration
3. **Permission errors:** Check file permissions in data directories

### Development Installation

For development, install additional tools:

```bash
pip install pytest black flake8 sphinx
```

### Virtual Environment (Recommended)

Create an isolated environment:

```bash
python -m venv rascas_env
source rascas_env/bin/activate  # On Windows: rascas_env\Scripts\activate
pip install -r requirements.txt
pip install -e .
```

## Next Steps

After installation:

1. Review the main [README.md](README.md) for usage examples
2. Check the [configuration guide](config/README.md) 
3. Run the example workflow: `python examples/example_workflow.py`
4. Process your simulation data using the provided scripts

## Support

For installation issues, please:

1. Check the troubleshooting section above
2. Verify all prerequisites are met
3. Open an issue on GitHub with detailed error messages