"""
Ramses Data Processor for Rascas Analysis

This module provides functionality to process Ramses simulation output
for analysis with Rascas tools.
"""

import numpy as np
import os
from pathlib import Path
from typing import Dict, List, Optional, Union
import logging

logger = logging.getLogger(__name__)


class RamsesDataProcessor:
    """
    Main processor for Ramses simulation output data.
    
    This class handles reading, processing and preparing Ramses output
    for Rascas analysis workflows.
    """
    
    def __init__(self, data_path: Union[str, Path], config: Optional[Dict] = None):
        """
        Initialize the Ramses data processor.
        
        Args:
            data_path: Path to Ramses output directory
            config: Optional configuration dictionary
        """
        self.data_path = Path(data_path)
        self.config = config or {}
        self.snapshots = []
        self._validate_data_path()
        
    def _validate_data_path(self):
        """Validate that the data path exists and contains Ramses output."""
        if not self.data_path.exists():
            raise FileNotFoundError(f"Data path does not exist: {self.data_path}")
            
        # Look for typical Ramses output files
        ramses_files = list(self.data_path.glob("output_*"))
        if not ramses_files:
            logger.warning(f"No Ramses output files found in {self.data_path}")
            
    def discover_snapshots(self) -> List[str]:
        """
        Discover available snapshots in the data directory.
        
        Returns:
            List of snapshot identifiers
        """
        snapshots = []
        
        # Look for Ramses output directories (output_00001, output_00002, etc.)
        for output_dir in sorted(self.data_path.glob("output_*")):
            if output_dir.is_dir():
                snapshot_id = output_dir.name.split("_")[-1]
                snapshots.append(snapshot_id)
                
        self.snapshots = snapshots
        logger.info(f"Found {len(snapshots)} snapshots: {snapshots}")
        return snapshots
        
    def load_snapshot(self, snapshot_id: str) -> Dict:
        """
        Load a specific snapshot for analysis.
        
        Args:
            snapshot_id: The snapshot identifier (e.g., "00001")
            
        Returns:
            Dictionary containing snapshot data
        """
        output_dir = self.data_path / f"output_{snapshot_id.zfill(5)}"
        
        if not output_dir.exists():
            raise FileNotFoundError(f"Snapshot directory not found: {output_dir}")
            
        # Load AMR and particle data
        snapshot_data = {
            'snapshot_id': snapshot_id,
            'path': output_dir,
            'amr_data': self._load_amr_data(output_dir),
            'particle_data': self._load_particle_data(output_dir),
            'info': self._load_info_file(output_dir)
        }
        
        return snapshot_data
        
    def _load_amr_data(self, output_dir: Path) -> Dict:
        """Load AMR (Adaptive Mesh Refinement) data."""
        # This would typically use specialized Ramses readers
        # For now, return placeholder structure
        amr_files = list(output_dir.glob("amr_*.out*"))
        return {
            'files': [str(f) for f in amr_files],
            'loaded': False,
            'grid_structure': None
        }
        
    def _load_particle_data(self, output_dir: Path) -> Dict:
        """Load particle data from Ramses output."""
        part_files = list(output_dir.glob("part_*.out*"))
        return {
            'files': [str(f) for f in part_files],
            'loaded': False,
            'particles': None
        }
        
    def _load_info_file(self, output_dir: Path) -> Dict:
        """Load the info file containing simulation parameters."""
        info_files = list(output_dir.glob("info_*.txt"))
        if not info_files:
            return {}
            
        info_file = info_files[0]
        info_data = {}
        
        try:
            with open(info_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if '=' in line:
                        key, value = line.split('=', 1)
                        info_data[key.strip()] = value.strip()
        except Exception as e:
            logger.warning(f"Could not parse info file {info_file}: {e}")
            
        return info_data
        
    def prepare_rascas_input(self, snapshot_data: Dict, output_path: Union[str, Path]) -> str:
        """
        Prepare Ramses data for Rascas analysis.
        
        Args:
            snapshot_data: Snapshot data from load_snapshot()
            output_path: Path where to save Rascas input files
            
        Returns:
            Path to the prepared Rascas input file
        """
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Create Rascas configuration for this snapshot
        rascas_config_file = output_path / f"rascas_config_{snapshot_data['snapshot_id']}.cfg"
        
        config_content = self._generate_rascas_config(snapshot_data)
        
        with open(rascas_config_file, 'w') as f:
            f.write(config_content)
            
        logger.info(f"Created Rascas configuration: {rascas_config_file}")
        return str(rascas_config_file)
        
    def _generate_rascas_config(self, snapshot_data: Dict) -> str:
        """Generate Rascas configuration content."""
        config_lines = [
            "# Rascas configuration for Ramses output",
            f"# Generated for snapshot: {snapshot_data['snapshot_id']}",
            "",
            "[input]",
            f"data_path = {snapshot_data['path']}",
            "format = ramses",
            "",
            "[analysis]",
            "compute_profiles = true",
            "compute_spectra = true", 
            "compute_statistics = true",
            "",
            "[output]",
            "save_plots = true",
            "save_data = true",
            ""
        ]
        
        return '\n'.join(config_lines)