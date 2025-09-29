"""
Gadget Data Processor for Rascas Analysis

This module provides functionality to process Gadget simulation output
for analysis with Rascas tools.
"""

import numpy as np
import os
from pathlib import Path
from typing import Dict, List, Optional, Union
import logging
import struct

logger = logging.getLogger(__name__)


class GadgetDataProcessor:
    """
    Main processor for Gadget simulation output data.
    
    This class handles reading, processing and preparing Gadget output
    for Rascas analysis workflows.
    """
    
    def __init__(self, data_path: Union[str, Path], config: Optional[Dict] = None):
        """
        Initialize the Gadget data processor.
        
        Args:
            data_path: Path to Gadget output directory
            config: Optional configuration dictionary
        """
        self.data_path = Path(data_path)
        self.config = config or {}
        self.snapshots = []
        self._validate_data_path()
        
    def _validate_data_path(self):
        """Validate that the data path exists and contains Gadget output."""
        if not self.data_path.exists():
            raise FileNotFoundError(f"Data path does not exist: {self.data_path}")
            
        # Look for typical Gadget output files
        gadget_files = list(self.data_path.glob("snapshot_*"))
        if not gadget_files:
            logger.warning(f"No Gadget snapshot files found in {self.data_path}")
            
    def discover_snapshots(self) -> List[str]:
        """
        Discover available snapshots in the data directory.
        
        Returns:
            List of snapshot identifiers
        """
        snapshots = []
        
        # Look for Gadget snapshot files (snapshot_000, snapshot_001, etc.)
        for snapshot_file in sorted(self.data_path.glob("snapshot_*")):
            if snapshot_file.is_file():
                snapshot_id = snapshot_file.name.split("_")[-1]
                snapshots.append(snapshot_id)
                
        self.snapshots = snapshots
        logger.info(f"Found {len(snapshots)} snapshots: {snapshots}")
        return snapshots
        
    def load_snapshot(self, snapshot_id: str) -> Dict:
        """
        Load a specific snapshot for analysis.
        
        Args:
            snapshot_id: The snapshot identifier (e.g., "000")
            
        Returns:
            Dictionary containing snapshot data
        """
        snapshot_file = self.data_path / f"snapshot_{snapshot_id.zfill(3)}"
        
        if not snapshot_file.exists():
            raise FileNotFoundError(f"Snapshot file not found: {snapshot_file}")
            
        # Load Gadget snapshot data
        snapshot_data = {
            'snapshot_id': snapshot_id,
            'path': snapshot_file,
            'header': self._load_header(snapshot_file),
            'particles': self._load_particles(snapshot_file)
        }
        
        return snapshot_data
        
    def _load_header(self, snapshot_file: Path) -> Dict:
        """Load Gadget snapshot header."""
        header_data = {}
        
        try:
            with open(snapshot_file, 'rb') as f:
                # Read Gadget header (simplified)
                # Real implementation would parse the full Gadget header format
                record_size = struct.unpack('I', f.read(4))[0]
                
                # Basic header fields (this is a simplified version)
                npart = struct.unpack('6I', f.read(24))  # Number of particles per type
                mass = struct.unpack('6d', f.read(48))   # Mass of particle types
                time = struct.unpack('d', f.read(8))[0]  # Time
                
                header_data = {
                    'npart': npart,
                    'mass': mass,
                    'time': time,
                    'record_size': record_size
                }
                
        except Exception as e:
            logger.warning(f"Could not parse Gadget header {snapshot_file}: {e}")
            header_data = {'error': str(e)}
            
        return header_data
        
    def _load_particles(self, snapshot_file: Path) -> Dict:
        """Load particle data from Gadget snapshot."""
        # This would typically use specialized Gadget readers like pynbody
        # For now, return placeholder structure
        return {
            'loaded': False,
            'positions': None,
            'velocities': None,
            'masses': None,
            'ids': None,
            'file': str(snapshot_file)
        }
        
    def prepare_rascas_input(self, snapshot_data: Dict, output_path: Union[str, Path]) -> str:
        """
        Prepare Gadget data for Rascas analysis.
        
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
            "# Rascas configuration for Gadget output",
            f"# Generated for snapshot: {snapshot_data['snapshot_id']}",
            "",
            "[input]",
            f"data_path = {snapshot_data['path']}",
            "format = gadget",
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
        
    def convert_to_hdf5(self, snapshot_data: Dict, output_file: Union[str, Path]) -> str:
        """
        Convert Gadget snapshot to HDF5 format for easier analysis.
        
        Args:
            snapshot_data: Snapshot data from load_snapshot()
            output_file: Output HDF5 file path
            
        Returns:
            Path to the created HDF5 file
        """
        try:
            import h5py
        except ImportError:
            raise ImportError("h5py is required for HDF5 conversion")
            
        output_file = Path(output_file)
        
        # This would implement the actual conversion
        # For now, create a placeholder HDF5 file
        with h5py.File(output_file, 'w') as h5f:
            h5f.attrs['snapshot_id'] = snapshot_data['snapshot_id']
            h5f.attrs['original_file'] = str(snapshot_data['path'])
            h5f.attrs['conversion_time'] = np.datetime64('now').astype(str)
            
            # Placeholder groups
            header_grp = h5f.create_group('Header')
            particles_grp = h5f.create_group('PartType0')  # Gas particles
            
            # Add header data
            if 'header' in snapshot_data:
                for key, value in snapshot_data['header'].items():
                    if isinstance(value, (int, float, str)):
                        header_grp.attrs[key] = value
                        
        logger.info(f"Created HDF5 file: {output_file}")
        return str(output_file)