"""
Plotting utilities for Rascas analysis
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)


class PlotManager:
    """Manager for creating analysis plots."""
    
    def __init__(self, config: Optional[Dict] = None):
        self.config = config or {}
        self._setup_style()
        
    def _setup_style(self):
        """Setup matplotlib style."""
        style = self.config.get('style', 'seaborn')
        try:
            plt.style.use(style)
        except Exception as e:
            logger.warning(f"Could not set style '{style}': {e}")
            
    def plot_density_profile(self, radius: np.ndarray, density: np.ndarray,
                           output_file: Optional[str] = None) -> str:
        """Plot density profile."""
        fig, ax = plt.subplots(figsize=self.config.get('figure_size', [10, 8]))
        
        ax.loglog(radius, density)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Density')
        ax.set_title('Density Profile')
        ax.grid(True, alpha=0.3)
        
        if output_file:
            plt.savefig(output_file, dpi=self.config.get('figure_dpi', 300))
            plt.close()
            return output_file
        else:
            plt.show()
            return ""