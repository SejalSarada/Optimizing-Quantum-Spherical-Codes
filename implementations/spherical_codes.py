"""
Spherical Code Optimization with Persistent Homology Constraints

This module implements the core functionality for optimizing high-dimensional
spherical codes using QAOA with topological constraints.

Author: Sejal Sarada
Institution: BITS Pilani, Goa Campus
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Union
from dataclasses import dataclass
import logging
from scipy.spatial.distance import pdist, squareform
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import AerSimulator
from qiskit.primitives import Estimator
from qiskit.algorithms.optimizers import COBYLA, SPSA

from .qaoa_optimizer import QAOAOptimizer
from .homology_constraints import PersistentHomologyConstraint

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class OptimizationResult:
    """Container for optimization results."""
    selected_points: np.ndarray
    min_distance: float
    betti_0_count: int
    betti_1_count: int
    logical_fidelity: float
    cost_evolution: List[float]
    optimal_parameters: Dict[str, np.ndarray]
    computation_time: float


class SphericalCodeOptimizer:
    """
    Main class for optimizing spherical codes with persistent homology constraints.
    
    This class implements the complete pipeline for:
    1. Generating candidate points on high-dimensional spheres
    2. Formulating the optimization problem
    3. Running QAOA with topological constraints
    4. Evaluating the resulting codes
    """
    
    def __init__(
        self,
        dimension: int,
        num_codewords: int,
        num_candidates: int = 20,
        constraint_weights: Optional[Dict[str, float]] = None,
        random_seed: Optional[int] = None
    ):
        """
        Initialize the spherical code optimizer.
        
        Args:
            dimension: Dimension of the sphere (d for S^{d-1})
            num_codewords: Number of codewords to select (N)
            num_candidates: Number of candidate points to generate (M)
            constraint_weights: Weights for persistent homology constraints
            random_seed: Random seed for reproducibility
        """
        self.dimension = dimension
        self.num_codewords = num_codewords
        self.num_candidates = num_candidates
        self.random_seed = random_seed
        
        # Default constraint weights
        if constraint_weights is None:
            constraint_weights = {'betti_0': 100.0, 'betti_1': 50.0}
        self.constraint_weights = constraint_weights
        
        # Initialize random number generator
        if random_seed is not None:
            np.random.seed(random_seed)
            
        # Generate candidate points
        self.candidate_points = self._generate_candidate_points()
        self.distance_matrix = self._compute_distance_matrix()
        
        # Initialize components
        self.qaoa_optimizer = None
        self.homology_constraint = PersistentHomologyConstraint()
        
        logger.info(f"Initialized SphericalCodeOptimizer: d={dimension}, N={num_codewords}, M={num_candidates}")
    
    def _generate_candidate_points(self) -> np.ndarray:
        """
        Generate candidate points uniformly distributed on the unit sphere.
        
        Returns:
            Array of shape (num_candidates, dimension) containing candidate points
        """
        # Generate points from standard Gaussian distribution
        points = np.random.randn(self.num_candidates, self.dimension)
        
        # Normalize to unit sphere
        norms = np.linalg.norm(points, axis=1, keepdims=True)
        points = points / norms
        
        return points
    
    def _compute_distance_matrix(self) -> np.ndarray:
        """
        Compute pairwise squared Euclidean distances between candidate points.
        
        Returns:
            Distance matrix of shape (num_candidates, num_candidates)
        """
        distances = pdist(self.candidate_points, metric='euclidean')
        return squareform(distances ** 2)  # Squared distances
    
    def _construct_cost_hamiltonian(self, include_homology: bool = True) -> SparsePauliOp:
        """
        Construct the cost Hamiltonian for QAOA optimization.
        
        Args:
            include_homology: Whether to include persistent homology constraints
            
        Returns:
            SparsePauliOp representing the cost Hamiltonian
        """
        # Distance-based cost terms
        pauli_list = []
        coeffs = []
        
        # Add pairwise distance terms: -∑_{i,j} D_{ij} x_i x_j
        for i in range(self.num_candidates):
            for j in range(i + 1, self.num_candidates):
                # Create Pauli string for x_i * x_j term
                pauli_str = ['I'] * self.num_candidates
                pauli_str[i] = 'Z'
                pauli_str[j] = 'Z'
                
                pauli_list.append(''.join(pauli_str))
                coeffs.append(-0.25 * self.distance_matrix[i, j])  # Factor of 0.25 for ZZ -> x_i x_j
        
        # Add constraint for fixed number of codewords: α(∑_i x_i - N)^2
        constraint_weight = 10.0
        
        # Linear terms: -2αN ∑_i x_i
        for i in range(self.num_candidates):
            pauli_str = ['I'] * self.num_candidates
            pauli_str[i] = 'Z'
            
            pauli_list.append(''.join(pauli_str))
            coeffs.append(constraint_weight * self.num_codewords * 0.5)  # Factor of 0.5 for Z -> x_i
        
        # Quadratic terms: α ∑_{i,j} x_i x_j
        for i in range(self.num_candidates):
            for j in range(i + 1, self.num_candidates):
                pauli_str = ['I'] * self.num_candidates
                pauli_str[i] = 'Z'
                pauli_str[j] = 'Z'
                
                pauli_list.append(''.join(pauli_str))
                coeffs.append(constraint_weight * 0.25)
        
        # Constant terms
        constant = constraint_weight * self.num_codewords ** 2
        
        return SparsePauliOp.from_list(list(zip(pauli_list, coeffs))) + constant
    
    def _construct_mixer_hamiltonian(self) -> SparsePauliOp:
        """
        Construct the mixer Hamiltonian for QAOA.
        
        Returns:
            SparsePauliOp representing the mixer Hamiltonian
        """
        pauli_list = []
        coeffs = []
        
        # Standard X mixer: ∑_i X_i
        for i in range(self.num_candidates):
            pauli_str = ['I'] * self.num_candidates
            pauli_str[i] = 'X'
            
            pauli_list.append(''.join(pauli_str))
            coeffs.append(1.0)
        
        return SparsePauliOp.from_list(list(zip(pauli_list, coeffs)))
    
    def _evaluate_solution(self, bitstring: str) -> Tuple[float, Dict[str, float]]:
        """
        Evaluate a candidate solution (bitstring) and compute metrics.
        
        Args:
            bitstring: Binary string representing selected points
            
        Returns:
            Tuple of (cost_value, metrics_dict)
        """
        selected_indices = [i for i, bit in enumerate(bitstring) if bit == '1']
        
        if len(selected_indices) != self.num_codewords:
            # Penalize invalid solutions
            return float('inf'), {}
        
        selected_points = self.candidate_points[selected_indices]
        
        # Compute minimum distance
        if len(selected_points) > 1:
            distances = pdist(selected_points, metric='euclidean')
            min_distance = np.min(distances)
        else:
            min_distance = 0.0
        
        # Compute persistent homology
        betti_0, betti_1 = self.homology_constraint.compute_betti_numbers(selected_points)
        
        # Compute cost including homology penalty
        distance_cost = -np.sum(self.distance_matrix[np.ix_(selected_indices, selected_indices)])
        homology_penalty = (
            self.constraint_weights['betti_0'] * (betti_0 - 1) ** 2 +
            self.constraint_weights['betti_1'] * betti_1 ** 2
        )
        
        total_cost = distance_cost + homology_penalty
        
        metrics = {
            'min_distance': min_distance,
            'betti_0': betti_0,
            'betti_1': betti_1,
            'distance_cost': distance_cost,
            'homology_penalty': homology_penalty,
            'total_cost': total_cost
        }
        
        return total_cost, metrics
    
    def _estimate_logical_fidelity(self, selected_points: np.ndarray, noise_level: float = 1e-3) -> float:
        """
        Estimate logical fidelity under depolarizing noise.
        
        Args:
            selected_points: Selected codeword points
            noise_level: Depolarizing noise probability
            
        Returns:
            Estimated logical fidelity
        """
        # Simplified model: fidelity decreases with smaller minimum distance
        if len(selected_points) > 1:
            distances = pdist(selected_points, metric='euclidean')
            min_distance = np.min(distances)
            
            # Heuristic: better separated codes have higher fidelity
            base_fidelity = 1 - noise_level * self.dimension
            distance_factor = min_distance / 2.0  # Normalize by maximum possible distance
            
            return base_fidelity * (1 + distance_factor)
        
        return 0.5  # Random guess for single point
    
    def optimize_with_constraints(
        self,
        qaoa_depth: int = 3,
        optimizer: str = 'COBYLA',
        max_iterations: int = 100,
        include_homology: bool = True
    ) -> OptimizationResult:
        """
        Run the complete optimization with persistent homology constraints.
        
        Args:
            qaoa_depth: Depth of QAOA ansatz (p parameter)
            optimizer: Classical optimizer to use
            max_iterations: Maximum number of optimization iterations
            include_homology: Whether to include homology constraints
            
        Returns:
            OptimizationResult containing optimization results
        """
        import time
        start_time = time.time()
        
        logger.info(f"Starting optimization with p={qaoa_depth}, optimizer={optimizer}")
        
        # Initialize QAOA optimizer
        cost_hamiltonian = self._construct_cost_hamiltonian(include_homology)
        mixer_hamiltonian = self._construct_mixer_hamiltonian()
        
        self.qaoa_optimizer = QAOAOptimizer(
            cost_hamiltonian=cost_hamiltonian,
            mixer_hamiltonian=mixer_hamiltonian,
            num_qubits=self.num_candidates,
            qaoa_depth=qaoa_depth
        )
        
        # Set up classical optimizer
        if optimizer == 'COBYLA':
            classical_optimizer = COBYLA(maxiter=max_iterations)
        elif optimizer == 'SPSA':
            classical_optimizer = SPSA(maxiter=max_iterations)
        else:
            raise ValueError(f"Unsupported optimizer: {optimizer}")
        
        # Run optimization
        result = self.qaoa_optimizer.optimize(
            optimizer=classical_optimizer,
            evaluation_callback=self._evaluation_callback
        )
        
        # Extract best solution
        best_bitstring = result.optimal_bitstring
        best_cost, best_metrics = self._evaluate_solution(best_bitstring)
        
        # Get selected points
        selected_indices = [i for i, bit in enumerate(best_bitstring) if bit == '1']
        selected_points = self.candidate_points[selected_indices]
        
        # Estimate logical fidelity
        logical_fidelity = self._estimate_logical_fidelity(selected_points)
        
        # Create result object
        optimization_result = OptimizationResult(
            selected_points=selected_points,
            min_distance=best_metrics.get('min_distance', 0.0),
            betti_0_count=best_metrics.get('betti_0', 0),
            betti_1_count=best_metrics.get('betti_1', 0),
            logical_fidelity=logical_fidelity,
            cost_evolution=result.cost_evolution,
            optimal_parameters=result.optimal_parameters,
            computation_time=time.time() - start_time
        )
        
        logger.info(f"Optimization completed in {optimization_result.computation_time:.2f}s")
        logger.info(f"Best min_distance: {optimization_result.min_distance:.4f}")
        logger.info(f"Betti numbers: β₀={optimization_result.betti_0_count}, β₁={optimization_result.betti_1_count}")
        
        return optimization_result
    
    def _evaluation_callback(self, cost: float, parameters: np.ndarray, iteration: int):
        """Callback function for monitoring optimization progress."""
        if iteration % 10 == 0:
            logger.info(f"Iteration {iteration}: cost = {cost:.6f}")
    
    def compare_methods(
        self,
        qaoa_depth: int = 3,
        num_random_trials: int = 10
    ) -> Dict[str, OptimizationResult]:
        """
        Compare different optimization methods.
        
        Args:
            qaoa_depth: QAOA depth for quantum methods
            num_random_trials: Number of random baseline trials
            
        Returns:
            Dictionary mapping method names to optimization results
        """
        results = {}
        
        # Random baseline
        logger.info("Running random baseline...")
        best_random_result = None
        best_random_distance = 0.0
        
        for trial in range(num_random_trials):
            # Randomly select codewords
            selected_indices = np.random.choice(
                self.num_candidates, 
                size=self.num_codewords, 
                replace=False
            )
            selected_points = self.candidate_points[selected_indices]
            
            # Evaluate
            if len(selected_points) > 1:
                distances = pdist(selected_points, metric='euclidean')
                min_distance = np.min(distances)
            else:
                min_distance = 0.0
            
            if min_distance > best_random_distance:
                best_random_distance = min_distance
                betti_0, betti_1 = self.homology_constraint.compute_betti_numbers(selected_points)
                logical_fidelity = self._estimate_logical_fidelity(selected_points)
                
                best_random_result = OptimizationResult(
                    selected_points=selected_points,
                    min_distance=min_distance,
                    betti_0_count=betti_0,
                    betti_1_count=betti_1,
                    logical_fidelity=logical_fidelity,
                    cost_evolution=[],
                    optimal_parameters={},
                    computation_time=0.0
                )
        
        results['random_baseline'] = best_random_result
        
        # QAOA unconstrained
        logger.info("Running QAOA unconstrained...")
        results['qaoa_unconstrained'] = self.optimize_with_constraints(
            qaoa_depth=qaoa_depth,
            include_homology=False
        )
        
        # QAOA with homology constraints
        logger.info("Running QAOA with homology constraints...")
        results['qaoa_constrained'] = self.optimize_with_constraints(
            qaoa_depth=qaoa_depth,
            include_homology=True
        )
        
        return results
    
    @classmethod
    def from_config(cls, config: Dict) -> 'SphericalCodeOptimizer':
        """
        Create optimizer from configuration dictionary.
        
        Args:
            config: Configuration parameters
            
        Returns:
            Configured SphericalCodeOptimizer instance
        """
        return cls(
            dimension=config['dimension'],
            num_codewords=config['num_codewords'],
            num_candidates=config.get('num_candidates', 20),
            constraint_weights=config.get('constraint_weights'),
            random_seed=config.get('random_seed')
        )
    
    def save_results(self, result: OptimizationResult, filename: str):
        """
        Save optimization results to file.
        
        Args:
            result: Optimization result to save
            filename: Output filename
        """
        import json
        
        # Convert result to serializable format
        result_dict = {
            'selected_points': result.selected_points.tolist(),
            'min_distance': result.min_distance,
            'betti_0_count': result.betti_0_count,
            'betti_1_count': result.betti_1_count,
            'logical_fidelity': result.logical_fidelity,
            'cost_evolution': result.cost_evolution,
            'computation_time': result.computation_time,
            'configuration': {
                'dimension': self.dimension,
                'num_codewords': self.num_codewords,
                'num_candidates': self.num_candidates,
                'constraint_weights': self.constraint_weights
            }
        }
        
        with open(filename, 'w') as f:
            json.dump(result_dict, f, indent=2)
        
        logger.info(f"Results saved to {filename}")


def main():
    """Example usage of the SphericalCodeOptimizer."""
    # Configuration for a simple 3D tetrahedral code
    config = {
        'dimension': 3,
        'num_codewords': 4,
        'num_candidates': 20,
        'constraint_weights': {'betti_0': 100.0, 'betti_1': 50.0},
        'random_seed': 42
    }
    
    # Initialize optimizer
    optimizer = SphericalCodeOptimizer.from_config(config)
    
    # Run comparison
    results = optimizer.compare_methods(qaoa_depth=3)
    
    # Print results
    print("\nComparison Results:")
    print("-" * 60)
    for method, result in results.items():
        print(f"{method:20s}: min_dist={result.min_distance:.4f}, "
              f"β₁={result.betti_1_count}, fidelity={result.logical_fidelity:.4f}")


if __name__ == "__main__":
    main()
