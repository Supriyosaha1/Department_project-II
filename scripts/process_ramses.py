#!/usr/bin/env python3
"""
Process Ramses simulation output for Rascas analysis

This script processes Ramses output files and prepares them for analysis
with Rascas tools.
"""

import argparse
import logging
from pathlib import Path
import sys

# Add the parent directory to the path so we can import our package
sys.path.insert(0, str(Path(__file__).parent.parent))

from rascas_analysis.ramses import RamsesDataProcessor
from rascas_analysis.utils import ConfigManager, RascasRunner


def setup_logging(verbose: bool = False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


def main():
    """Main processing function."""
    parser = argparse.ArgumentParser(description='Process Ramses output for Rascas analysis')
    parser.add_argument('data_path', help='Path to Ramses output directory')
    parser.add_argument('--output', '-o', default='./rascas_output', 
                       help='Output directory for Rascas configurations')
    parser.add_argument('--config', '-c', help='Configuration file')
    parser.add_argument('--snapshots', nargs='+', help='Specific snapshots to process')
    parser.add_argument('--run-rascas', action='store_true', 
                       help='Run Rascas analysis after preparing configs')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    logger = logging.getLogger(__name__)
    logger.info("Starting Ramses data processing")
    
    # Load configuration
    config_manager = ConfigManager(args.config) if args.config else ConfigManager()
    
    # Initialize Ramses processor
    try:
        processor = RamsesDataProcessor(args.data_path, config_manager.get_section('ramses'))
        logger.info(f"Initialized Ramses processor for: {args.data_path}")
    except Exception as e:
        logger.error(f"Failed to initialize Ramses processor: {e}")
        return 1
    
    # Discover snapshots
    available_snapshots = processor.discover_snapshots()
    if not available_snapshots:
        logger.error("No snapshots found in data directory")
        return 1
    
    # Determine which snapshots to process
    if args.snapshots:
        snapshots_to_process = [s for s in args.snapshots if s in available_snapshots]
        if not snapshots_to_process:
            logger.error(f"None of the specified snapshots found: {args.snapshots}")
            return 1
    else:
        snapshots_to_process = available_snapshots
    
    logger.info(f"Processing {len(snapshots_to_process)} snapshots: {snapshots_to_process}")
    
    # Process each snapshot
    config_files = []
    output_path = Path(args.output)
    
    for snapshot_id in snapshots_to_process:
        try:
            logger.info(f"Processing snapshot {snapshot_id}")
            
            # Load snapshot data
            snapshot_data = processor.load_snapshot(snapshot_id)
            
            # Prepare Rascas input
            config_file = processor.prepare_rascas_input(
                snapshot_data, 
                output_path / f"snapshot_{snapshot_id}"
            )
            config_files.append(config_file)
            
            logger.info(f"Created Rascas configuration for snapshot {snapshot_id}")
            
        except Exception as e:
            logger.error(f"Error processing snapshot {snapshot_id}: {e}")
            continue
    
    if not config_files:
        logger.error("No configurations created successfully")
        return 1
    
    logger.info(f"Successfully created {len(config_files)} Rascas configurations")
    
    # Optionally run Rascas analysis
    if args.run_rascas:
        logger.info("Running Rascas analysis...")
        runner = RascasRunner(config_manager=config_manager)
        results = runner.run_batch_analysis(config_files, output_path / "analysis_results")
        
        successful = sum(1 for r in results if r['success'])
        logger.info(f"Rascas analysis completed: {successful}/{len(results)} successful")
    
    logger.info("Processing completed successfully")
    return 0


if __name__ == '__main__':
    sys.exit(main())