# Configuration Guide

This directory contains configuration files and templates for the Rascas analysis framework.

## Configuration Files

### `default.json`
Main configuration file containing all default settings. This file uses JSON format and includes:

- **general**: Global settings (debug mode, parallelization, etc.)
- **ramses**: Ramses-specific settings (unit systems, data loading options)
- **gadget**: Gadget-specific settings (particle types, formats)
- **rascas**: Rascas tool configuration (executable path, output formats)
- **analysis**: Analysis parameters (binning, computations to perform)
- **plotting**: Visualization settings (figure formats, styles, DPI)

### `rascas_template.cfg`
Template configuration file for Rascas analysis. This file is used as a base when generating Rascas configuration files for specific snapshots.

## Configuration Sections

### General Settings
```json
{
  "general": {
    "debug": false,
    "verbose": true,
    "parallel": true,
    "num_threads": 4
  }
}
```

### Ramses Configuration
```json
{
  "ramses": {
    "default_unit_system": "cgs",
    "load_particles": true,
    "load_amr": true,
    "max_level": null
  }
}
```

### Gadget Configuration
```json
{
  "gadget": {
    "default_unit_system": "gadget",
    "load_gas": true,
    "load_dark_matter": true,
    "load_stars": true
  }
}
```

### Rascas Tool Settings
```json
{
  "rascas": {
    "executable": "rascas",
    "config_template": "default.cfg",
    "output_format": "hdf5"
  }
}
```

## Using Configuration Files

### From Command Line
Specify a configuration file using the `--config` option:

```bash
python scripts/process_ramses.py /path/to/data --config config/my_config.json
```

### From Python Code
```python
from rascas_analysis.utils import ConfigManager

# Load specific configuration
config = ConfigManager("config/default.json")

# Access configuration sections
ramses_config = config.get_section('ramses')
analysis_settings = config.get_section('analysis')

# Get specific values
debug_mode = config.get('debug', section='general', default=False)
```

## Creating Custom Configurations

### Method 1: Copy and Modify
```bash
cp config/default.json config/my_analysis.json
# Edit my_analysis.json with your settings
```

### Method 2: Create Programmatically
```python
from rascas_analysis.utils import ConfigManager

config = ConfigManager()
config.set('executable', '/custom/path/to/rascas', 'rascas')
config.set('num_threads', 8, 'general')
config.save_config('config/custom.json')
```

## Environment Variables

Some settings can be overridden using environment variables:

- `RASCAS_EXECUTABLE`: Path to Rascas executable
- `RASCAS_THREADS`: Number of threads to use
- `RASCAS_DEBUG`: Enable debug mode (true/false)

## Validation

The configuration system includes validation for:

- File paths existence
- Numerical ranges
- Required parameters
- Format compatibility

## Examples

See the `examples/` directory for complete workflow examples using different configurations.

## Troubleshooting

### Common Configuration Issues

1. **Invalid JSON syntax**: Use a JSON validator or Python's json module
2. **Missing executable paths**: Verify Rascas installation and update paths
3. **Incorrect data types**: Ensure numbers are not quoted in JSON
4. **Path separators**: Use forward slashes or escape backslashes in paths

### Debugging Configuration

Enable debug mode to see detailed configuration loading:

```json
{
  "general": {
    "debug": true,
    "verbose": true
  }
}
```