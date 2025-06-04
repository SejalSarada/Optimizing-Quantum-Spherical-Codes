"""
QAOA Optimizer for Spherical Codes with Persistent Homology Constraints

This module implements the Quantum Approximate Optimization Algorithm (QAOA)
for optimizing spherical code configurations with topological constraints.

Author: Sejal Sarada
BITS Pilani, K.K. Birla Goa Campus
"""

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.algorithms.optimizers import COBYLA
from qiskit.providers.aer import AerSimulator
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import PauliEvolutionGate
from qiskit import transpile, execute
from typing import List, Tuple, Dict, Optional
import warnings
warnings.filterwarnings('ignore')

from .homology_constraints import PersistentHomologyCalculator
from .spherical_codes import SphericalCodeGenerator


class QAOASphericalOptimizer:
    """
    QAOA optimizer for spherical codes with persistent homology constraints.
    
    This class implements the core QAOA algorithm with custom cost and mixing
    Hamiltonians designed for spherical code optimization.
    """
    
    def __init__(self, 
                 dimension: int, 
                 num_codewords: int,
                 num_candidates: int,
                 p_layers: int = 3,
                 shots: int = 1024,
                 lambda_0: float = 10.0,
                 lambda_1: float = 10.0):
        """
        Initialize the QAOA optimizer.
        
        Args:
            dimension: Dimension of the sphere (d for S^{d-1})
            num_codewords: Number of codewords to select (N)
            num_candidates: Number of candidate points (M)
            p_layers: Number of QAOA layers
            shots: Number of measurement shots
            lambda_0: Weight for Betti_0 constraint
            lambda_1: Weight for Betti_1 constraint
        """
        self.dimension = dimension
        self.num_codewords = num_codewords
        self.num_candidates = num_candidates
        self.p_layers = p_layers
        self.shots = shots
        self.lambda_0 = lambda_0
        self.lambda_1 = lambda_1
        
        # Initialize components
        self.sphere_generator = SphericalCodeGenerator(dimension)
        self.homology_calculator = PersistentHomologyCalculator()
        self.backend = AerSimulator()
        
        # Generate candidate points
        self.candidate_points = self.sphere_generator.generate_random_points(num_candidates)
        self.distance_matrix = self._compute_distance_matrix()
        
        # Optimization history
        self.optimization_history = []
        
    def _compute_distance_matrix(self) -> np.ndarray:
        """Compute pairwise squared Euclidean distances between candidate points."""
        M = len(self.candidate_points)
        D = np.zeros((M, M))
        
        for i in range(M):
            for j in range(i+1, M):
                dist_sq = np.sum((self.candidate_points[i] - self.candidate_points[j])**2)
                D[i, j] = D[j, i] = dist_sq
                
        return D
    
    def _create_cost_hamiltonian(self) -> SparsePauliOp:
        """
        Create the cost Hamiltonian H_C = -sum_{i,j} D_{ij} * x_i * x_j
        """
        pauli_list = []
        coeffs = []
        
        # Distance terms
        for i in range(self.num_candidates):
            for j in range(i+1, self.num_candidates):
                if self.distance_matrix[i, j] > 0:
                    # Create Pauli string for Z_i * Z_j interaction
                    pauli_str = ['I'] * self.num_candidates
                    pauli_str[i] = 'Z'
                    pauli_str[j] = 'Z'
                    
                    pauli_list.append(''.join(pauli_str))
                    coeffs.append(-self.distance_matrix[i, j] / 4)  # Factor of 1/4 for Z->x conversion
        
        # Constraint penalty: (sum_i x_i - N)^2
        constraint_weight = 5.0
        
        # Linear terms: -2N * sum_i x_i
        for i in range(self.num_candidates):
            pauli_str = ['I'] * self.num_candidates
            pauli_str[i] = 'Z'
            pauli_list.append(''.join(pauli_str))
            coeffs.append(constraint_weight * self.num_codewords / 2)
        
        # Quadratic terms: sum_i sum_j x_i x_j
        for i in range(self.num_candidates):
            for j in range(i+1, self.num_candidates):
                pauli_str = ['I'] * self.num_candidates
                pauli_str[i] = 'Z'
                pauli_str[j] = 'Z'
                pauli_list.append(''.join(pauli_str))
                coeffs.append(-constraint_weight / 4)
        
        # Constant term (handled separately)
        constant = constraint_weight * (self.num_codewords**2 + self.num_candidates / 2)
        
        return SparsePauliOp(pauli_list, coeffs)
    
    def _create_mixing_hamiltonian(self) -> SparsePauliOp:
        """Create the mixing Hamiltonian H_B = sum_i X_i"""
        pauli_list = []
        coeffs = []
        
        for i in range(self.num_candidates):
            pauli_str = ['I'] * self.num_candidates
            pauli_str[i] = 'X'
            pauli_list.append(''.join(pauli_str))
            coeffs.append(1.0)
            
        return SparsePauliOp(pauli_list, coeffs)
    
    def _create_qaoa_circuit(self, params: np.ndarray) -> QuantumCircuit:
        """
        Create the QAOA circuit with given parameters.
        
        Args:
            params: Array of 2*p parameters [gamma_1, beta_1, ..., gamma_p, beta_p]
        """
        qreg = QuantumRegister(self.num_candidates, 'q')
        creg = ClassicalRegister(self.num_candidates, 'c')
        circuit = QuantumCircuit(qreg, creg)
        
        # Initialize in superposition
        circuit.h(range(self.num_candidates))
        
        # Apply QAOA layers
        cost_hamiltonian = self._create_cost_hamiltonian()
        mixing_hamiltonian = self._create_mixing_hamiltonian()
        
        for layer in range(self.p_layers):
            gamma = params[2*layer]
            beta = params[2*layer + 1]
            
            # Cost unitary
            cost_evolution = PauliEvolutionGate(cost_hamiltonian, time=gamma)
            circuit.append(cost_evolution, range(self.num_candidates))
            
            # Mixing unitary
            mixing_evolution = PauliEvolutionGate(mixing_hamiltonian, time=beta)
            circuit.append(mixing_evolution, range(self.num_candidates))
        
        # Measurement
        circuit.measure(range(self.num_candidates), range(self.num_candidates))
        
        return circuit
    
    def _evaluate_bitstring(self, bitstring: str) -> Tuple[float, Dict]:
        """
        Evaluate the cost of a given bitstring solution.
        
        Args:
            bitstring: Binary string representing point selection
            
        Returns:
            Tuple of (cost, metrics_dict)
        """
        # Convert bitstring to selection indices
        selected_indices = [i for i, bit in enumerate(bitstring) if bit == '1']
        
        if len(selected_indices) != self.num_codewords:
            # Penalty for wrong number of codewords
            return float('inf'), {}
        
        # Get selected points
        selected_points = [self.candidate_points[i] for i in selected_indices]
        
        # Calculate minimum distance
        min_distance = float('inf')
        for i in range(len(selected_points)):
            for j in range(i+1, len(selected_points)):
                dist = np.linalg.norm(selected_points[i] - selected_points[j])
                min_distance = min(min_distance, dist)
        
        # Calculate persistent homology metrics
        betti_0, betti_1 = self.homology_calculator.compute_betti_numbers(
            selected_points, threshold=0.5
        )
        
        # Compute cost
        distance_cost = -min_distance  # Negative because we want to maximize
        homology_penalty = self.lambda_0 * (betti_0 - 1)**2 + self.lambda_1 * betti_1**2
        
        total_cost = distance_cost + homology_penalty
        
        metrics = {
            'min_distance': min_distance,
            'betti_0': betti_0,
            'betti_1': betti_1,
            'distance_cost': distance_cost,
            'homology_penalty': homology_penalty,
            'selected_indices': selected_indices
        }
        
        return total_cost, metrics
    
    def _objective_function(self, params: np.ndarray) -> float:
        """
        Objective function for classical optimization.
        
        Args:
            params: QAOA parameters
            
        Returns:
            Expected cost value
        """
        # Create and execute circuit
        circuit = self._create_qaoa_circuit(params)
        compiled_circuit = transpile(circuit, self.backend)
        job = execute(compiled_circuit, self.backend, shots=self.shots)
        result = job.result()
        counts = result.get_counts()
        
        # Calculate expected cost
        total_cost = 0.0
        total_counts = 0
        
        for bitstring, count in counts.items():
            cost, metrics = self._evaluate_bitstring(bitstring)
            if cost != float('inf'):
                total_cost += cost * count
                total_counts += count
        
        if total_counts == 0:
            return float('inf')
        
        expected_cost = total_cost / total_counts
        
        # Store optimization history
        self.optimization_history.append({
            'params': params.copy(),
            'cost': expected_cost,
            'counts': counts
        })
        
        return expected_cost
    
    def optimize(self, 
                 initial_params: Optional[np.ndarray] = None,
                 maxiter: int = 100) -> Dict:
        """
        Run the QAOA optimization.
        
        Args:
            initial_params: Initial parameter values
            maxiter: Maximum optimization iterations
            
        Returns:
            Dictionary containing optimization results
        """
        # Initialize parameters
        if initial_params is None:
            initial_params = np.random.uniform(0, 2*np.pi, 2*self.p_layers)
        
        # Set up optimizer
        optimizer = COBYLA(maxiter=maxiter)
        
        # Run optimization
        print(f"Starting QAOA optimization with {self.p_layers} layers...")
        result = optimizer.minimize(
            fun=self._objective_function,
            x0=initial_params
        )
        
        # Get best solution
        best_params = result.x
        best_cost = result.fun
        
        # Analyze best solution
        circuit = self._create_qaoa_circuit(best_params)
        compiled_circuit = transpile(circuit, self.backend)
        job = execute(compiled_circuit, self.backend, shots=self.shots)
        final_result = job.result()
        final_counts = final_result.get_counts()
        
        # Find the most frequent valid bitstring
        best_bitstring = None
        best_metrics = None
        max_count = 0
        
        for bitstring, count in final_counts.items():
            cost, metrics = self._evaluate_bitstring(bitstring)
            if cost != float('inf') and count > max_count:
                max_count = count
                best_bitstring = bitstring
                best_metrics = metrics
        
        return {
            'optimal_params': best_params,
            'optimal_cost': best_cost,
            'best_bitstring': best_bitstring,
            'best_metrics': best_metrics,
            'final_counts': final_counts,
            'optimization_history': self.optimization_history,
            'convergence_info': {
                'success': result.success,
                'nfev': result.nfev,
                'message': result.message
            }
        }
    
    def get_selected_points(self, bitstring: str) -> np.ndarray:
        """Get the actual selected points from a bitstring."""
        selected_indices = [i for i, bit in enumerate(bitstring) if bit == '1']
        return np.array([self.candidate_points[i] for i in selected_indices])
