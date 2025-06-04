"""
Persistent Homology Constraints for Spherical Codes

This module implements persistent homology calculations using the Vietoris-Rips complex
to enforce topological constraints on spherical code configurations.

Author: Sejal Sarada
BITS Pilani, K.K. Birla Goa Campus
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
import warnings
warnings.filterwarnings('ignore')

try:
    import gudhi as gd
    GUDHI_AVAILABLE = True
except ImportError:
    GUDHI_AVAILABLE = False
    warnings.warn("GUDHI not available. Using simplified homology calculations.")


class PersistentHomologyCalculator:
    """
    Calculator for persistent homology features of point clouds.
    
    This class computes Vietoris-Rips complexes and extracts persistent Betti numbers
    to characterize the topological structure of spherical code configurations.
    """
    
    def __init__(self, max_dimension: int = 2):
        """
        Initialize the persistent homology calculator.
        
        Args:
            max_dimension: Maximum dimension for homology computation
        """
        self.max_dimension = max_dimension
        self.persistence_diagrams = {}
        
    def compute_vietoris_rips_complex(self, 
                                    points: np.ndarray, 
                                    max_edge_length: float) -> Optional[object]:
        """
        Compute the Vietoris-Rips complex for given points.
        
        Args:
            points: Array of shape (n_points, dimension)
            max_edge_length: Maximum edge length for complex construction
            
        Returns:
            Vietoris-Rips complex object (if GUDHI available)
        """
        if not GUDHI_AVAILABLE:
            return None
            
        # Create Vietoris-Rips complex
        rips_complex = gd.RipsComplex(points=points, max_edge_length=max_edge_length)
        simplex_tree = rips_complex.create_simplex_tree(max_dimension=self.max_dimension)
        
        return simplex_tree
    
    def compute_persistence(self, simplex_tree) -> List[Tuple]:
        """
        Compute persistence pairs from a simplex tree.
        
        Args:
            simplex_tree: GUDHI simplex tree object
            
        Returns:
            List of persistence pairs (dimension, (birth, death))
        """
        if not GUDHI_AVAILABLE or simplex_tree is None:
            return []
            
        # Compute persistence
        persistence = simplex_tree.persistence()
        return persistence
    
    def extract_betti_numbers(self, 
                            persistence: List[Tuple], 
                            threshold: float) -> Tuple[int, int]:
        """
        Extract Betti numbers at a given threshold.
        
        Args:
            persistence: List of persistence pairs
            threshold: Threshold value for persistence
            
        Returns:
            Tuple of (betti_0, betti_1)
        """
        if not persistence:
            return 1, 0  # Default for disconnected case
            
        betti_0 = 0  # Connected components
        betti_1 = 0  # 1-dimensional holes (loops)
        
        for dim, (birth, death) in persistence:
            # Check if feature persists at threshold
            if birth <= threshold and (death > threshold or death == float('inf')):
                if dim == 0:
                    betti_0 += 1
                elif dim == 1:
                    betti_1 += 1
                    
        return betti_0, betti_1
    
    def compute_betti_numbers(self, 
                            points: np.ndarray, 
                            threshold: float = 0.5) -> Tuple[int, int]:
        """
        Main method to compute Betti numbers for a point cloud.
        
        Args:
            points: Array of shape (n_points, dimension) representing point cloud
            threshold: Threshold for persistence computation
            
        Returns:
            Tuple of (betti_0, betti_1)
        """
        if len(points) < 2:
            return 1, 0
        
        if GUDHI_AVAILABLE:
            # Use GUDHI for accurate computation
            try:
                # Compute maximum distance for complex construction
                max_dist = self._compute_max_distance(points)
                max_edge_length = min(threshold * 2, max_dist)
                
                # Build Vietoris-Rips complex
                simplex_tree = self.compute_vietoris_rips_complex(points, max_edge_length)
                
                if simplex_tree is None:
                    return self._simplified_betti_computation(points, threshold)
                
                # Compute persistence
                persistence = self.compute_persistence(simplex_tree)
                
                # Extract Betti numbers
                betti_0, betti_1 = self.extract_betti_numbers(persistence, threshold)
                
                # Store persistence diagram for analysis
                self.persistence_diagrams[id(points)] = persistence
                
                return betti_0, betti_1
                
            except Exception as e:
                warnings.warn(f"GUDHI computation failed: {e}. Using simplified method.")
                return self._simplified_betti_computation(points, threshold)
        else:
            return self._simplified_betti_computation(points, threshold)
    
    def _simplified_betti_computation(self, 
                                    points: np.ndarray, 
                                    threshold: float) -> Tuple[int, int]:
        """
        Simplified Betti number computation without GUDHI.
        
        This method provides approximate Betti numbers based on distance analysis.
        """
        n_points = len(points)
        
        # Compute pairwise distances
        distances = np.zeros((n_points, n_points))
        for i in range(n_points):
            for j in range(i+1, n_points):
                dist = np.linalg.norm(points[i] - points[j])
                distances[i, j] = distances[j, i] = dist
        
        # Estimate connected components (Betti_0)
        # Use Union-Find to count connected components at threshold
        parent = list(range(n_points))
        
        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]
        
        def union(x, y):
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py
        
        # Connect points within threshold distance
        for i in range(n_points):
            for j in range(i+1, n_points):
                if distances[i, j] <= threshold:
                    union(i, j)
        
        # Count connected components
        components = len(set(find(i) for i in range(n_points)))
        betti_0 = components
        
        # Estimate loops (Betti_1) - simplified heuristic
        # Count potential triangles that could indicate loops
        betti_1 = 0
        if n_points >= 3:
            triangle_threshold = threshold * 1.5
            triangles = 0
            
            for i in range(n_points):
                for j in range(i+1, n_points):
                    for k in range(j+1, n_points):
                        if (distances[i, j] <= triangle_threshold and 
                            distances[j, k] <= triangle_threshold and 
                            distances[i, k] <= triangle_threshold):
                            triangles += 1
            
            # Heuristic: estimate loops based on triangle density
            if triangles > n_points:
                betti_1 = max(0, triangles - n_points + 1)
        
        return betti_0, betti_1
    
    def _compute_max_distance(self, points: np.ndarray) -> float:
        """Compute maximum pairwise distance in point cloud."""
        max_dist = 0.0
        n_points = len(points)
        
        for i in range(n_points):
            for j in range(i+1, n_points):
                dist = np.linalg.norm(points[i] - points[j])
                max_dist = max(max_dist, dist)
                
        return max_dist
    
    def analyze_persistence_diagram(self, 
                                  points: np.ndarray) -> Dict[str, any]:
        """
        Analyze the persistence diagram for given points.
        
        Args:
            points: Point cloud array
            
        Returns:
            Dictionary with persistence analysis
        """
        # Compute Betti numbers at multiple thresholds
        thresholds = np.linspace(0.1, 2.0, 20)
        betti_evolution = {'threshold': [], 'betti_0': [], 'betti_1': []}
        
        for threshold in thresholds:
            b0, b1 = self.compute_betti_numbers(points, threshold)
            betti_evolution['threshold'].append(threshold)
            betti_evolution['betti_0'].append(b0)
            betti_evolution['betti_1'].append(b1)
        
        # Find persistent features
        persistent_features = self._find_persistent_features(betti_evolution)
        
        return {
            'betti_evolution': betti_evolution,
            'persistent_features': persistent_features,
            'summary': {
                'max_components': max(betti_evolution['betti_0']),
                'max_loops': max(betti_evolution['betti_1']),
                'stable_threshold': self._find_stable_threshold(betti_evolution)
            }
        }
    
    def _find_persistent_features(self, betti_evolution: Dict) -> Dict:
        """Find features that persist across multiple thresholds."""
        thresholds = betti_evolution['threshold']
        b0_values = betti_evolution['betti_0']
        b1_values = betti_evolution['betti_1']
        
        # Find threshold ranges where Betti numbers are stable
        stable_ranges = {'betti_0': [], 'betti_1': []}
        
        # For Betti_0
        current_value = b0_values[0]
        start_idx = 0
        
        for i in range(1, len(b0_values)):
            if b0_values[i] != current_value:
                if i - start_idx > 2:  # At least 3 consecutive points
                    stable_ranges['betti_0'].append({
                        'value': current_value,
                        'start_threshold': thresholds[start_idx],
                        'end_threshold': thresholds[i-1],
                        'persistence': thresholds[i-1] - thresholds[start_idx]
                    })
                current_value = b0_values[i]
                start_idx = i
        
        # Similar for Betti_1
        current_value = b1_values[0]
        start_idx = 0
        
        for i in range(1, len(b1_values)):
            if b1_values[i] != current_value:
                if i - start_idx > 2:
                    stable_ranges['betti_1'].append({
                        'value': current_value,
                        'start_threshold': thresholds[start_idx],
                        'end_threshold': thresholds[i-1],
                        'persistence': thresholds[i-1] - thresholds[start_idx]
                    })
                current_value = b1_values[i]
                start_idx = i
        
        return stable_ranges
    
    def _find_stable_threshold(self, betti_evolution: Dict) -> float:
        """Find a threshold where both Betti numbers are stable."""
        # Look for a threshold where Betti_0 = 1 and Betti_1 = 0
        for i, (b0, b1) in enumerate(zip(betti_evolution['betti_0'], 
                                        betti_evolution['betti_1'])):
            if b0 == 1 and b1 == 0:
                return betti_evolution['threshold'][i]
        
        # If not found, return middle threshold
        return betti_evolution['threshold'][len(betti_evolution['threshold'])//2]
    
    def visualize_persistence_diagram(self, points: np.ndarray, save_path: str = None):
        """
        Visualize the persistence diagram (requires matplotlib).
        
        Args:
            points: Point cloud array
            save_path: Optional path to save the plot
        """
        try:
            import matplotlib.pyplot as plt
            
            analysis = self.analyze_persistence_diagram(points)
            betti_evolution = analysis['betti_evolution']
            
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
            
            # Plot Betti number evolution
            ax1.plot(betti_evolution['threshold'], betti_evolution['betti_0'], 
                    'o-', label='Betti_0 (Components)', linewidth=2)
            ax1.plot(betti_evolution['threshold'], betti_evolution['betti_1'], 
                    's-', label='Betti_1 (Loops)', linewidth=2)
            ax1.set_xlabel('Threshold')
            ax1.set_ylabel('Betti Number')
            ax1.set_title('Betti Number Evolution')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            # Plot point cloud (for low dimensions)
            if points.shape[1] == 2:
                ax2.scatter(points[:, 0], points[:, 1], s=100, alpha=0.7)
                ax2.set_title('Point Cloud (2D projection)')
                ax2.set_aspect('equal')
            elif points.shape[1] == 3:
                ax2 = fig.add_subplot(122, projection='3d')
                ax2.scatter(points[:, 0], points[:, 1], points[:, 2], s=100, alpha=0.7)
                ax2.set_title('Point Cloud (3D)')
            else:
                ax2.text(0.5, 0.5, f'Point cloud in {points.shape[1]}D\n({points.shape[0]} points)',
                        ha='center', va='center', transform=ax2.transAxes, fontsize=12)
                ax2.set_title('High-dimensional Point Cloud')
            
            plt.tight_layout()
            
            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
            plt.show()
            
        except ImportError:
            print("Matplotlib not available for visualization.")
