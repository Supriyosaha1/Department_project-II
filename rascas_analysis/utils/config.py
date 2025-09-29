"""
Configuration management for Rascas analysis workflows
"""

import configparser
import json
from pathlib import Path
from typing import Dict, Any, Optional, Union
import logging

logger = logging.getLogger(__name__)


class ConfigManager:
    """
    Manages configuration for Rascas analysis workflows.
    
    Supports both INI-style config files and JSON configurations.
    """
    
    def __init__(self, config_file: Optional[Union[str, Path]] = None):
        """
        Initialize configuration manager.
        
        Args:
            config_file: Optional path to configuration file
        """
        self.config_file = Path(config_file) if config_file else None
        self.config = {}
        self.defaults = self._get_defaults()
        
        if self.config_file and self.config_file.exists():
            self.load_config()
        else:
            self.config = self.defaults.copy()
            
    def _get_defaults(self) -> Dict[str, Any]:
        """Get default configuration values."""
        return {
            'general': {
                'debug': False,
                'verbose': True,
                'parallel': True,
                'num_threads': 4
            },
            'ramses': {
                'default_unit_system': 'cgs',
                'load_particles': True,
                'load_amr': True,
                'max_level': None
            },
            'gadget': {
                'default_unit_system': 'gadget',
                'load_gas': True,
                'load_dark_matter': True,
                'load_stars': True
            },
            'rascas': {
                'executable': 'rascas',
                'config_template': 'default.cfg',
                'output_format': 'hdf5'
            },
            'analysis': {
                'compute_profiles': True,
                'compute_spectra': True,
                'compute_statistics': True,
                'radial_bins': 50,
                'mass_bins': 50
            },
            'plotting': {
                'figure_format': 'png',
                'figure_dpi': 300,
                'figure_size': [10, 8],
                'style': 'seaborn'
            }
        }
        
    def load_config(self, config_file: Optional[Union[str, Path]] = None):
        """
        Load configuration from file.
        
        Args:
            config_file: Path to configuration file (uses self.config_file if None)
        """
        if config_file:
            self.config_file = Path(config_file)
            
        if not self.config_file or not self.config_file.exists():
            logger.warning(f"Configuration file not found: {self.config_file}")
            return
            
        file_extension = self.config_file.suffix.lower()
        
        try:
            if file_extension == '.json':
                self._load_json_config()
            elif file_extension in ['.cfg', '.ini', '.conf']:
                self._load_ini_config()
            else:
                logger.warning(f"Unsupported config file format: {file_extension}")
                
        except Exception as e:
            logger.error(f"Error loading config file {self.config_file}: {e}")
            
    def _load_json_config(self):
        """Load JSON configuration file."""
        with open(self.config_file, 'r') as f:
            loaded_config = json.load(f)
            
        # Merge with defaults
        self.config = self._merge_configs(self.defaults, loaded_config)
        
    def _load_ini_config(self):
        """Load INI-style configuration file."""
        parser = configparser.ConfigParser()
        parser.read(self.config_file)
        
        loaded_config = {}
        for section_name in parser.sections():
            loaded_config[section_name] = {}
            for key, value in parser.items(section_name):
                # Try to convert to appropriate type
                loaded_config[section_name][key] = self._convert_value(value)
                
        # Merge with defaults
        self.config = self._merge_configs(self.defaults, loaded_config)
        
    def _convert_value(self, value: str) -> Any:
        """Convert string value to appropriate type."""
        # Handle boolean values
        if value.lower() in ['true', 'yes', '1']:
            return True
        elif value.lower() in ['false', 'no', '0']:
            return False
        elif value.lower() in ['none', 'null']:
            return None
            
        # Try to convert to number
        try:
            if '.' in value:
                return float(value)
            else:
                return int(value)
        except ValueError:
            pass
            
        # Try to parse as list
        if value.startswith('[') and value.endswith(']'):
            try:
                return json.loads(value)
            except json.JSONDecodeError:
                pass
                
        # Return as string
        return value
        
    def _merge_configs(self, default: Dict, loaded: Dict) -> Dict:
        """Merge loaded configuration with defaults."""
        merged = default.copy()
        
        for key, value in loaded.items():
            if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
                merged[key] = self._merge_configs(merged[key], value)
            else:
                merged[key] = value
                
        return merged
        
    def save_config(self, config_file: Optional[Union[str, Path]] = None, format: str = 'json'):
        """
        Save current configuration to file.
        
        Args:
            config_file: Output configuration file path
            format: File format ('json' or 'ini')
        """
        if config_file:
            output_file = Path(config_file)
        else:
            output_file = self.config_file
            
        if not output_file:
            raise ValueError("No output file specified")
            
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        if format == 'json':
            with open(output_file, 'w') as f:
                json.dump(self.config, f, indent=2)
        elif format == 'ini':
            parser = configparser.ConfigParser()
            for section_name, section_data in self.config.items():
                parser.add_section(section_name)
                for key, value in section_data.items():
                    parser.set(section_name, key, str(value))
                    
            with open(output_file, 'w') as f:
                parser.write(f)
        else:
            raise ValueError(f"Unsupported format: {format}")
            
        logger.info(f"Configuration saved to: {output_file}")
        
    def get(self, key: str, section: Optional[str] = None, default: Any = None) -> Any:
        """
        Get configuration value.
        
        Args:
            key: Configuration key
            section: Configuration section (if None, searches all sections)
            default: Default value if key not found
            
        Returns:
            Configuration value
        """
        if section:
            return self.config.get(section, {}).get(key, default)
        else:
            # Search all sections
            for section_data in self.config.values():
                if isinstance(section_data, dict) and key in section_data:
                    return section_data[key]
            return default
            
    def set(self, key: str, value: Any, section: str):
        """
        Set configuration value.
        
        Args:
            key: Configuration key
            value: Configuration value
            section: Configuration section
        """
        if section not in self.config:
            self.config[section] = {}
        self.config[section][key] = value
        
    def get_section(self, section: str) -> Dict[str, Any]:
        """
        Get entire configuration section.
        
        Args:
            section: Section name
            
        Returns:
            Section configuration dictionary
        """
        return self.config.get(section, {})