# Optimizing-Quantum-Spherical-Codes
Hybrid quantum-classical framework to optimize high-dimensional spherical codes using QAOA with persistent homology constraints. Demonstrates topologically regular quantum code configurations with improved fidelity under realistic noise models.

> **Author**: Sejal Sarada  
> **Advisors**: Dr. Sumit Kale (IBM Quantum India), Dr. Himadri Mukherjee (BITS Pilani Goa Campus)  
> **Thesis Date**: May 2025

---

## Overview  
**A Novel Framework for Quantum Error Correction Code Design using Quantum Spherical Codes with Persistent Homology Constraints**

This repository contains the implementation and results from the research thesis "Optimizing High Dimensional Spherical Codes with Persistent Homology Constraints using QAOA". We present a groundbreaking approach to quantum error-correcting code design that combines the Quantum Approximate Optimization Algorithm (QAOA) with topological data analysis through persistent homology constraints. Our method treats quantum codewords as points on high-dimensional spheres and optimizes their configurations to maximize error correction capabilities while enforcing desirable topological properties. This work was conducted as part of a Bachelor's thesis at BITS Pilani, Goa Campus, carried out under supervision of Dr. Sumit Kale.

---

## Key Innovation
This work represents the first integration of persistent homology constraints into quantum spherical code optimization, creating a novel intersection of:

- Quantum Error Correction (QEC)  
- Topological Data Analysis (TDA) 
- Quantum Optimization Algorithms 
- High-Dimensional Quantum Spherical code Geometry

---

## Research Contributions

### 1. Mathematical Framework
- **Topological Constraints**: Integration of Vietoris-Rips complexes and Betti number analysis into spherical code optimization 
- **Hard Constraint Enforcement**: Custom QAOA mixing Hamiltonians for maintaining fixed code sizes  
- **Multi-objective Optimization**: Simultaneous optimization of Euclidean distance and topological properties    

### 2. Algorithmic Innovation
- **PH-QAOA**: Novel QAOA variant incorporating persistent homology penalties
- **Constraint Hamiltonian**: H'_C = H_C + H_PH where H_PH penalizes undesirable topological features
- **Scalable Implementation**: Efficient computation for high-dimensional spherical codes (d=3-6)

### 3. Performance Achievements
- **140% improvement** in minimum chord distance over random placement
- **Complete elimination** of spurious topological holes (Betti₁ = 0)
- **15-20% enhancement** in logical fidelity under noise
- **Near-optimal configurations** approaching theoretical bounds 

---

## Technical Highlights

### Spherical Code Optimization
- **Cost Hamiltonian**: H_C = -∑_{i,j} D_ij x_i x_j  
- **Topological Penalty**: H_PH = λ₀(b₀ - 1)² + λ₁b₁²  
- **Combined Objective**: H'_C = H_C + H_PH

### Persistent Homology Integration
- **Betti₀ constraint**: Ensures single connected component (b₀ = 1)  
- **Betti₁ constraint**: Eliminates persistent loops (b₁ = 0)  
- **Multi-scale analysis**: Vietoris-Rips complex filtration  
- **Topological persistence**: Robust feature detection across scales  

### QAOA Implementation
- **Depth optimization**: Systematic evaluation of p = 1–4 layers  
- **Custom mixers**: Hard constraint preservation during evolution  
- **Hybrid optimization**: Classical parameter tuning with COBYLA  
- **Quantum simulation**: Qiskit Aer with 1024 shots per evaluation  

## Experimental Results
The following have been plotted, displayed and elaborated on in the report OptimizingQuantumSphericalCodes_ThesisReport.pdf

### Code Performance Comparison

| Configuration         | Method              | Min. Distance | Betti₁ Count | Logical Fidelity |
|-----------------------|---------------------|----------------|----------------|-------------------|
| d=3, N=4              | Random Baseline     | 0.66           | 0.38           | 0.20              |
| d=3, N=4              | QAOA (unconstrained)| 1.60           | 0.08           | 0.59              |
| d=3, N=4              | PH-QAOA             | 1.58           | 0.00           | 0.73              |
| d=4, N=5              | Random Baseline     | 0.52           | 0.72           | 0.17              |
| d=4, N=5              | QAOA (unconstrained)| 1.25           | 0.36           | 0.58              |
| d=4, N=5              | PH-QAOA             | 1.20           | 0.00           | 0.65              |

### Key Findings
- Topological regularization consistently eliminates spurious homological features  
- Distance optimization approaches theoretical maxima (e.g., ~1.63 for tetrahedral codes)  
- Error correction enhancement through improved geometric uniformity  
- Scalability demonstrated across multiple dimensions and code sizes  

## Implementation Details

### Software Stack
- **Quantum Framework**: IBM Qiskit v0.28  
- **Simulation Backend**: Qiskit Aer (statevector and noisy simulation)  
- **Topological Analysis**: GUDHI 3.7.1, Ripser++ 1.1.0  
- **Classical Optimization**: COBYLA with 100 function evaluations  
- **Hardware Models**: FakeMelbourne, FakeAlmaden for realistic noise  

### Parameter Configuration
- **Dimensions**: d = 3–6 (spheres S^{d-1})  
- **Code sizes**: N = 3–6 codewords  
- **Candidate sets**: M = 20–30 points per sphere  
- **QAOA depth**: p = 1–4 layers  
- **Noise model**: Depolarizing errors (~10⁻³ rate)  

## Applications and Impact

### Quantum Error Correction
- **Bosonic codes**: Direct application to multimode quantum systems  
- **High-dimensional QEC**: Scalable approach for large Hilbert spaces  
- **Fault-tolerant computing**: Foundation for robust quantum computation  

### Algorithmic Advances
- **Constrained QAOA**: Framework for topological optimization problems  
- **TDA-guided optimization**: Novel application of persistent homology  
- **Quantum-classical hybrid**: Efficient integration of classical TDA with quantum optimization  

## Future Directions

### Near-term Development
- Hardware validation on NISQ devices (IBM, Rigetti platforms)  
- Error mitigation integration for practical quantum implementations  
- Advanced topological features (Stiefel-Whitney classes, zigzag persistence)  
- Automated feature selection for problem-specific optimization  

### Long-term Vision
- Fully quantum TDA: Integration with quantum persistent homology algorithms  
- Dynamical codes: Extension to time-varying error correction schemes  
- Machine learning integration: Data-driven code optimization  
- Theoretical bounds: Rigorous performance guarantees for topologically-constrained codes

---

## Repository Structure

universal-dataset-encoder/ <br />
├── README.md <br />
├── thesis-report/ <br />
│ └── OptimizingQuantumSphericalCodes_ThesisReport.pdf <br />
└── implementations/ <br />
├── qaoa_spherical.py      # QAOA optimization routines <br />
└── persistent_homology.py # TDA constraint computation <br />
└── spherical_codes.py     # Core spherical code implementatio <br />

---

## Citation

If you use this work in your research, please cite:

bibtex <br />
@thesis{sarada2025universal, <br />
&nbsp; title={Optimizing High Dimensional Spherical Codes with Persistent Homology Constraints using QAOA}, <br />
&nbsp;author={Sarada, Sejal}, <br />
&nbsp;year={2025}, <br />
&nbsp;school={BITS Pilani, Goa Campus}, <br />
&nbsp;supervisor={Sumit Suresh Kale and Himadri Mukherjee},<br />
&nbsp;type={Bachelor's Thesis}<br />
}

---

## Contributors

### Research & Development  
- **Sejal Sarada** — BITS Pilani, Goa Campus ; [email: sejalsarada13@gmail.com]

### Supervision  
- **Dr. Sumit Suresh Kale** — Qiskit Software Developer and Quantum Researcher, IBM Quantum India ; [email: sumit.suresh.kale1@ibm.com]
- **Dr. Himadri Mukherjee** — Assistant Professor, BITS Pilani Goa Campus ; [email: himadrim@goa.bits-pilani.ac.in]

---

## Contact

For inquiries or collaboration opportunities:

- **Research inquiries**: Reach out through academic email [sejalsarada13@gmail.com] or institutional channels  
- **Technical issues**: [Open an issue](https://github.com/SejalSarada/Optimizing-Quantum-Spherical-Codes/issues)  
- **Collaborations**: Connect via [LinkedIn]([https://www.linkedin.com/in/sejalsarada/](https://www.linkedin.com/in/sejal-sarada-88ab96204/))

---

> This repository contains the complete implementation of the research conducted for the Bachelor's thesis "Optimizing High Dimensional Spherical Codes with Persistent Homology Constraints using QAOA" submitted to BITS Pilani in May 2025.
