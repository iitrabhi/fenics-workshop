# Introduction to FEniCS

This is a five-day course focused on solving partial differential equations (PDEs) using the FEniCS software package. The goal is to introduce the students to PDEs encountered in various engineering and science disciplines, such as solid mechanics, heat transfer, and mass transport. 

The course materials, including tutorials and exercises, were created as part of a five-day workshop at IIT Madras, in collaboration with Vanderbilt University, USA. These materials are presented in [Jupyter Notebooks](https://jupyter.org/), which allow you to see both the code and its explanations, as well as the results, all together.

The `tutorials` are comprehensive notebooks that demonstrate how to approach different types of problems using FEniCS. On the other hand, the `exercises` are meant to be interactive, and they encourage you to expand the notebooks by adding new functionalities. This way, you can develop your expertise in using FEniCS.

If you just want to view the tutorials and exercises without making any changes, you can access them on [on nbviewer](https://nbviewer.org/github/iitrabhi/iitm-fenics-course/blob/48a6a14f8f7c27f2a32cf1ea101e18934d254b01/README.ipynb) without installing FEniCS. However, if you want to edit and undertake the exercises, you will need to  [install FEniCS](install-instructions.ipynb). Additionally, you have the option to either clone the repository or download the code from the provided link.

## What is FEniCS

FEniCS is an acronym that stands for "Finite Element Computational Software." The inclusion of "ni" in the name is to create a balanced and appealing composition. The FEniCS software package was compiled at the University of Chicago, whose Phoenix mascot likely influenced the choice of the name.

FEniCS is a high-performance computing (HPC) capable tool that efficiently utilizes supercomputers and high-performance clusters to solve complex scientific problems. It supports parallel computing, JIT compilation, and integrates with PETSc and MPI for scalability and performance. Its HPC capabilities enable researchers to perform large-scale simulations and analyses effectively.

## Contents
- [Instructions for installing FEniCS](install-instructions.ipynb)
- Day 1
  - Tutorials
    - [Tutorial 1-1: Solving linear Poisson's equation](src/day-1/tutorials/1_linear_poisson.ipynb)
    - [Tutorial 1-2: Visualization](src/day-1/tutorials/2_visualization.ipynb)
  - Exercises
    - [Exercise 1-1: Built in meshes](src/day-1/exercises/1_built_in_mesh.ipynb)
    - [Exercise 1-2: Boundary conditions](src/day-1/exercises/2_boundary_conditions.ipynb)
    - [Exercise 1-3: Expressions](src/day-1/exercises/3_expressions.ipynb)
    - [Exercise 1-4: Spacially varying properties](src/day-1/exercises/4_spacially_varying_properties.ipynb)
    - [Exercise 1-5: Convergence](src/day-1/exercises/5_convergence.ipynb)
- Day 2
  - Tutorials
    - [Tutorial 2-1: Solving non-linear Poisson's equation using Picard Iteration](src/day-2/tutorials/1_non_linear_poisson_picard.ipynb)
    - [Tutorial 2-3: Solving non-linear Poisson's equation using Newton Iteration](src/day-2/tutorials/2_non_linear_poisson_newton.ipynb)
    - [Tutorial 2-3: Newton method with Manual Differentiation](src/day-2/tutorials/3_non_linear_poisson_newton_manual_diff.ipynb)
    - [Tutorial 2-4: Newton method with Automatic Differentiation](src/day-2/tutorials/4_non_linear_poisson_newton_auto_diff.ipynb)
  - Exercises
    - [Exercise 2-1: 2D and 3D Domains](src/day-2/exercises/1_2d_3d_domains.ipynb)
    - [Exercise 2-2: Boundary Conditions](src/day-2/exercises/2_boundary_conditions.ipynb)
    - [Exercise 2-3: Non-Linearity and Tolerances](src/day-2/exercises/3_non_linearity_tolerances.ipynb)
    - [Exercise 2-4: Manual vs Auto differentiation](src/day-2/exercises/4_manual_auto_differentiation.ipynb)
- Day 3
  - Tutorials
    - [Tutorial 3-1: Beam bending](src/day-3/tutorials/1_beam_bending.ipynb)
    - [Tutorial 3-2: Loads and boundary conditions](src/day-3/tutorials/2_load_and_boundary_conditions.ipynb)
    - [Tutorial 3-3: Solver design](src/day-3/tutorials/3_solver.ipynb)
  - Exercises
    - [Exercise 3-1: 3D Uniaxial](src/day-3/exercises/1_3d_uniaxial.ipynb)
    - [Exercise 3-2: Plane Stress v/s Plane strain](src/day-3/exercises/2_plane_stress_strain.ipynb)
    - [Exercise 3-3: Effect of solver design](src/day-3/exercises/3_solver_design.ipynb)
- Day 4
  - Tutorials
    - [Tutorial 4-1: Solving hyper elasticity in FEniCS](src/day-4/tutorials/1_hyper_elastic.ipynb)
    - [Tutorial 4-1: Pseudo time stepping in hyper elasticity](src/day-4/tutorials/2_load_displacement.ipynb)
  - Exercises
    - [Exercise 4-1: Mooney rivlin material model](src/day-4/exercises/1_mooney-rivlin.ipynb)
    - [Exercise 4-2: Comparison between elastic and hyperelastic solution](src/day-4/exercises/2_comparison.ipynb)
- Day 5
  - Tutorials
    - [Tutorial 5-1:Thermoelasticity](src/day-5/tutorials/1_thermoelasticity.ipynb)
    - [Tutorial 5-2:Multi-material domains](src/day-5/tutorials/2_multi_material.ipynb)
    - [Tutorial 5-3:Transient analysis](src/day-5/tutorials/3_transient_analysis.ipynb)
  - Exercise
    - [Exercise 5-1:Bi-metallic quasi static analysis](src/day-5/exercises/1_bi_metallic_quasi_static.ipynb)
    - [Exercise 5-2:Transient thermo-mechanical analysis](src/day-5/exercises/2_thermo_mechanical_transient.ipynb)



## Where to find help

When in doubt about any specific FEniCS command or implementation, there are several resources you can refer to for clarification and assistance:

1. **FEniCS Book:** The official FEniCS book is an invaluable resource. It provides comprehensive documentation, tutorials, examples, and detailed explanations of all FEniCS functionalities. You can find it at: [The FEniCS Book](https://fenicsproject.org/book/)
2. **FEniCS Q&A Forum:** The FEniCS Q&A forum is an active community where users and developers discuss issues, ask questions, and share knowledge about FEniCS. You can post your questions here and get help from experienced FEniCS users and developers. Access the forum at [FEniCS discourse](https://fenicsproject.discourse.group/)
3. **FEniCS Workshop Slack channel:** You can also join the FEniCS workshop community on Slack to connect with other users, developers, and experts in the field. Slack is an excellent platform for real-time discussions, sharing knowledge, and getting help with FEniCS-related topics. You can join the Slack community by using [this link](https://join.slack.com/t/fenicsworkshop/shared_invite/zt-1zqrdlvmr-LziJQ59NdEOBzn7YZf59hg). 
4. **FEniCS online documentation:** The online documentation since V1.3.0 is available at [this link](https://fenicsproject.org/olddocs/dolfin/). You can go into the python folder of any version to find the syntax and meaning of different commands. The latest version is not that descriptive, and thus I recommend checking out older documentation, as most of the commands are the same.

Note: Please go through the official community [page here](https://fenicsproject.org/community/).

## Books to read
1. Bleyer, Jeremy. “Numerical Tours of Continuum Mechanics Using FEniCS,” n.d., 105.
2. Langtangen, Hans Petter, and Anders Logg. “Solving PDEs in Python – The FEniCS Tutorial Volume I,” n.d., 153.
3. Langtangen, Hans Petter, and Kent-Andre Mardal. “Introduction to Numerical Methods for Variational Problems”. Vol. 21. Texts in Computational Science and Engineering. Cham: Springer International Publishing, 2019. https://doi.org/10.1007/978-3-030-23788-2.
4. Logg, Anders, Kent-Andre Mardal, and Garth Wells, eds. “Automated Solution of Differential Equations by the Finite Element Method”. Vol. 84. Lecture Notes in Computational Science and Engineering. Berlin, Heidelberg: Springer Berlin Heidelberg, 2012. https://doi.org/10.1007/978-3-642-23099-8.
