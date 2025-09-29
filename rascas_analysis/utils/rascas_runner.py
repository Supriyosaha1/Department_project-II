"""
Rascas runner utility for executing Rascas analysis workflows
"""

import subprocess
import os
from pathlib import Path
from typing import Dict, List, Optional, Union
import logging
import time

logger = logging.getLogger(__name__)


class RascasRunner:
    """
    Utility class for running Rascas analysis on processed simulation data.
    """
    
    def __init__(self, rascas_executable: str = "rascas", config_manager=None):
        """
        Initialize Rascas runner.
        
        Args:
            rascas_executable: Path to Rascas executable
            config_manager: Optional ConfigManager instance
        """
        self.rascas_executable = rascas_executable
        self.config_manager = config_manager
        self._validate_executable()
        
    def _validate_executable(self):
        """Validate that Rascas executable is available."""
        try:
            result = subprocess.run(
                [self.rascas_executable, "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                logger.info(f"Rascas executable found: {self.rascas_executable}")
                logger.info(f"Version info: {result.stdout.strip()}")
            else:
                logger.warning(f"Rascas executable returned error: {result.stderr}")
        except FileNotFoundError:
            logger.warning(f"Rascas executable not found: {self.rascas_executable}")
        except subprocess.TimeoutExpired:
            logger.warning("Rascas executable validation timed out")
        except Exception as e:
            logger.warning(f"Error validating Rascas executable: {e}")
            
    def run_analysis(self, config_file: Union[str, Path], 
                    output_dir: Optional[Union[str, Path]] = None,
                    verbose: bool = True) -> Dict:
        """
        Run Rascas analysis with given configuration.
        
        Args:
            config_file: Path to Rascas configuration file
            output_dir: Optional output directory (overrides config)
            verbose: Enable verbose output
            
        Returns:
            Dictionary containing run results and metadata
        """
        config_file = Path(config_file)
        if not config_file.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_file}")
            
        # Prepare command
        cmd = [self.rascas_executable]
        
        if verbose:
            cmd.append("--verbose")
            
        cmd.extend(["--config", str(config_file)])
        
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            cmd.extend(["--output", str(output_dir)])
            
        # Run analysis
        start_time = time.time()
        logger.info(f"Running Rascas analysis: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=config_file.parent
            )
            
            end_time = time.time()
            runtime = end_time - start_time
            
            analysis_result = {
                'success': result.returncode == 0,
                'runtime': runtime,
                'command': ' '.join(cmd),
                'config_file': str(config_file),
                'output_dir': str(output_dir) if output_dir else None,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'returncode': result.returncode
            }
            
            if result.returncode == 0:
                logger.info(f"Rascas analysis completed successfully in {runtime:.2f} seconds")
            else:
                logger.error(f"Rascas analysis failed with return code {result.returncode}")
                logger.error(f"Error output: {result.stderr}")
                
            return analysis_result
            
        except Exception as e:
            logger.error(f"Error running Rascas analysis: {e}")
            return {
                'success': False,
                'error': str(e),
                'command': ' '.join(cmd),
                'config_file': str(config_file)
            }
            
    def run_batch_analysis(self, config_files: List[Union[str, Path]],
                          output_base_dir: Optional[Union[str, Path]] = None,
                          parallel: bool = False) -> List[Dict]:
        """
        Run Rascas analysis on multiple configurations.
        
        Args:
            config_files: List of configuration file paths
            output_base_dir: Base directory for outputs
            parallel: Whether to run analyses in parallel (not implemented yet)
            
        Returns:
            List of analysis results
        """
        results = []
        
        for i, config_file in enumerate(config_files):
            logger.info(f"Running analysis {i+1}/{len(config_files)}: {config_file}")
            
            if output_base_dir:
                config_path = Path(config_file)
                output_dir = Path(output_base_dir) / f"analysis_{config_path.stem}"
            else:
                output_dir = None
                
            result = self.run_analysis(config_file, output_dir)
            results.append(result)
            
            if not result['success']:
                logger.warning(f"Analysis {i+1} failed, continuing with next...")
                
        successful = sum(1 for r in results if r['success'])
        logger.info(f"Batch analysis completed: {successful}/{len(results)} successful")
        
        return results
        
    def create_analysis_script(self, config_files: List[Union[str, Path]],
                              script_path: Union[str, Path],
                              output_base_dir: Optional[Union[str, Path]] = None) -> str:
        """
        Create a shell script for running multiple Rascas analyses.
        
        Args:
            config_files: List of configuration file paths
            script_path: Output script path
            output_base_dir: Base directory for outputs
            
        Returns:
            Path to created script
        """
        script_path = Path(script_path)
        script_lines = [
            "#!/bin/bash",
            "# Rascas batch analysis script",
            f"# Generated by RascasRunner",
            "",
            "set -e  # Exit on error",
            ""
        ]
        
        for i, config_file in enumerate(config_files):
            config_path = Path(config_file)
            
            if output_base_dir:
                output_dir = Path(output_base_dir) / f"analysis_{config_path.stem}"
                mkdir_cmd = f"mkdir -p {output_dir}"
                run_cmd = f"{self.rascas_executable} --config {config_path} --output {output_dir}"
            else:
                mkdir_cmd = ""
                run_cmd = f"{self.rascas_executable} --config {config_path}"
                
            script_lines.extend([
                f"echo \"Running analysis {i+1}/{len(config_files)}: {config_path}\"",
                mkdir_cmd if mkdir_cmd else "",
                run_cmd,
                "echo \"Analysis completed\"",
                ""
            ])
            
        script_content = '\n'.join(line for line in script_lines if line)
        
        with open(script_path, 'w') as f:
            f.write(script_content)
            
        # Make script executable
        script_path.chmod(0o755)
        
        logger.info(f"Created analysis script: {script_path}")
        return str(script_path)