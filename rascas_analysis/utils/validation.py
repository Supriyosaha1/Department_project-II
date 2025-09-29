"""
Data validation utilities
"""

from pathlib import Path
from typing import Dict, Any, List
import logging

logger = logging.getLogger(__name__)


class DataValidator:
    """Validator for simulation data."""
    
    @staticmethod
    def validate_ramses_output(data_path: Path) -> Dict[str, Any]:
        """Validate Ramses output directory."""
        validation_result = {
            'valid': True,
            'errors': [],
            'warnings': []
        }
        
        if not data_path.exists():
            validation_result['valid'] = False
            validation_result['errors'].append(f"Path does not exist: {data_path}")
            
        # Check for output directories
        output_dirs = list(data_path.glob("output_*"))
        if not output_dirs:
            validation_result['warnings'].append("No output directories found")
            
        return validation_result
        
    @staticmethod
    def validate_gadget_output(data_path: Path) -> Dict[str, Any]:
        """Validate Gadget output directory."""
        validation_result = {
            'valid': True,
            'errors': [],
            'warnings': []
        }
        
        if not data_path.exists():
            validation_result['valid'] = False
            validation_result['errors'].append(f"Path does not exist: {data_path}")
            
        # Check for snapshot files
        snapshot_files = list(data_path.glob("snapshot_*"))
        if not snapshot_files:
            validation_result['warnings'].append("No snapshot files found")
            
        return validation_result