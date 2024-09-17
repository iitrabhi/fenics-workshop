![](attachments/Pasted%20image%2020231219133930.png)
## Introduction to FEniCS

This is a five-day course focused on solving partial differential equations (PDEs) using the FEniCS software package. The goal is to introduce the students to PDEs encountered in various engineering and science disciplines, such as solid mechanics, heat transfer, and mass transport. 

The course materials, including tutorials and exercises, were created as part of a five-day workshop at IIT Madras, in collaboration with Vanderbilt University, USA. These materials are presented in [Jupyter Notebooks](https://jupyter.org/), which allow you to see both the code and its explanations, as well as the results, all together.

The `tutorials` are comprehensive notebooks that demonstrate how to approach different types of problems using FEniCS. On the other hand, the `exercises` are meant to be interactive, and they encourage you to expand the notebooks by adding new functionalities. This way, you can develop your expertise in using FEniCS.



## What is FEniCS?
FEniCS is a high-performance computing (HPC) capable tool that efficiently utilizes supercomputers and high-performance clusters to solve complex scientific problems. It supports parallel computing, JIT compilation, and integrates with PETSc and MPI for scalability and performance. Its HPC capabilities enable researchers to perform large-scale simulations and analyses effectively.

!!! note "FEniCS is an acronym"
	FEniCS is an acronym that stands for "Finite Element Computational Software." The inclusion of "ni" in the name is to create a balanced and appealing composition. The FEniCS software package was compiled at the University of Chicago, whose Phoenix mascot likely influenced the choice of the name.



## Contents 
 [Instructions for installing FEniCS](1_introduction/2_installation.md)


- Day 1
	- Tutorials
	    - [Tutorial 1-1: Solving linear Poisson's equation](2_notebooks/day-1/tutorials/1_linear_poisson.md)
	    - [Tutorial 1-2: Visualization](2_notebooks/day-1/tutorials/2_visualization.md)
	  - Exercises
	    - [Exercise 1-1: Built in meshes](2_notebooks/day-1/exercises/1_built_in_mesh.md)
	    - [Exercise 1-2: Boundary conditions](2_notebooks/day-1/exercises/2_boundary_conditions.md)
	    - [Exercise 1-3: Expressions](2_notebooks/day-1/exercises/3_expressions.md)
	    - [Exercise 1-4: Spacially varying properties](2_notebooks/day-1/exercises/4_spacially_varying_properties.md)
	    - [Exercise 1-5: Convergence](2_notebooks/day-1/exercises/5_convergence.md)
- Day 2
	- Tutorials
	    - [Tutorial 2-1: Solving non-linear Poisson's equation using Picard Iteration](2_notebooks/day-2/tutorials/1_non_linear_poisson_picard.md)
	    - [Tutorial 2-3: Solving non-linear Poisson's equation using Newton Iteration](2_notebooks/day-2/tutorials/2_non_linear_poisson_newton.md)
	    - [Tutorial 2-3: Newton method with Manual Differentiation](2_notebooks/day-2/tutorials/3_non_linear_poisson_newton_manual_diff.md)
	    - [Tutorial 2-4: Newton method with Automatic Differentiation](2_notebooks/day-2/tutorials/4_non_linear_poisson_newton_auto_diff.md)
	  - Exercises
	    - [Exercise 2-1: 2D and 3D Domains](2_notebooks/day-2/exercises/1_2d_3d_domains.md)
	    - [Exercise 2-2: Boundary Conditions](2_notebooks/day-2/exercises/2_boundary_conditions.md)
	    - [Exercise 2-3: Non-Linearity and Tolerances](2_notebooks/day-2/exercises/3_non_linearity_tolerances.md)
	    - [Exercise 2-4: Manual vs Auto differentiation](2_notebooks/day-2/exercises/4_manual_auto_differentiation.md)
- Day 3
	- Tutorials
	    - [Tutorial 3-1: Beam bending](2_notebooks/day-3/tutorials/1_beam_bending.md)
	    - [Tutorial 3-2: Loads and boundary conditions](2_notebooks/day-3/tutorials/2_load_and_boundary_conditions.md)
	    - [Tutorial 3-3: Solver design](2_notebooks/day-3/tutorials/3_solver.md)
	  - Exercises
	    - [Exercise 3-1: 3D Uniaxial](2_notebooks/day-3/exercises/1_3d_uniaxial.md)
	    - [Exercise 3-2: Plane Stress v/s Plane strain](2_notebooks/day-3/exercises/2_plane_stress_strain.md)
	    - [Exercise 3-3: Effect of solver design](2_notebooks/day-3/exercises/3_solver_design.md)
- Day 4
	- Tutorials
	    - [Tutorial 4-1: Solving hyper elasticity in FEniCS](2_notebooks/day-4/tutorials/1_hyper_elastic.md)
	    - [Tutorial 4-1: Pseudo time stepping in hyper elasticity](2_notebooks/day-4/tutorials/2_load_displacement.md)
	  - Exercises
	    - [Exercise 4-1: Mooney rivlin material model](2_notebooks/day-4/exercises/1_mooney-rivlin.md)
	    - [Exercise 4-2: Comparison between elastic and hyperelastic solution](2_notebooks/day-4/exercises/2_comparison.md)
- Day 5
	- Tutorials
	    - [Tutorial 5-1:Thermoelasticity](2_notebooks/day-5/tutorials/1_thermoelasticity.md)
	    - [Tutorial 5-2:Multi-material domains](2_notebooks/day-5/tutorials/2_multi_material.md)
	    - [Tutorial 5-3:Transient analysis](2_notebooks/day-5/tutorials/3_transient_analysis.md)
	  - Exercise
	    - [Exercise 5-1:Bi-metallic quasi static analysis](2_notebooks/day-5/exercises/1_bi_metallic_quasi_static.md)
	    - [Exercise 5-2:Transient thermo-mechanical analysis](2_notebooks/day-5/exercises/2_thermo_mechanical_transient.md) 

## Github Repo
- [FEniCS Workshop](https://github.com/iitrabhi/fenics-workshop)

## Authors
- [Abhinav Gupta](https://abhigupta.io)
- [Ravindra Duddu](https://engineering.vanderbilt.edu/bio/ravindra-duddu)


