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
    - [Tutorial 1-2: Visualization](src\day-1\tutorials\2_visualization.ipynb)
  - Exercises
    - [Exercise 1-1: Built in meshes](src\day-1\exercises\1_built_in_mesh.ipynb)
    - [Exercise 1-2: Boundary conditions](src\day-1\exercises\2_boundary_conditions.ipynb)
    - [Exercise 1-3: Expressions](src\day-1\exercises\3_expressions.ipynb)
    - [Exercise 1-4: Spacially varying properties](src\day-1\exercises\4_spacially_varying_properties.ipynb)
    - [Exercise 1-5: Convergence](src\day-1\exercises\5_convergence.ipynb)
- Day 2
  - Tutorials
    - [Tutorial 2-1: Solving non-linear Poisson's equation](src/day-2/tutorials/1_non_linear_poisson.ipynb)
    - [Tutorial 2-2: Auto differentiation](src\day-2\tutorials\2_auto-differentiation.ipynb)
    - [Tutorial 2-3: Solver design](src\day-2\tutorials\2_solver_design.ipynb)
  - Exercises
    - [Exercise 2-1: Built in meshes](src\day-1\exercises\1_built_in_mesh.ipynb)
    - [Exercise 2-2: Boundary conditions](src\day-1\exercises\2_boundary_conditions.ipynb)

## Where to find help

When in doubt about any specific FEniCS command or implementation, there are several resources you can refer to for clarification and assistance:

1. **FEniCS Book:** The official FEniCS book is an invaluable resource. It provides comprehensive documentation, tutorials, examples, and detailed explanations of all FEniCS functionalities. You can find it at: [The FEniCS Book](https://fenicsproject.org/book/)
2. **FEniCS Q&A Forum:** The FEniCS Q&A forum is an active community where users and developers discuss issues, ask questions, and share knowledge about FEniCS. You can post your questions here and get help from experienced FEniCS users and developers. Access the forum at [FEniCS discourse](https://fenicsproject.discourse.group/)
3. **FEniCS Workshop Slack channel:** You can also join the FEniCS workshop community on Slack to connect with other users, developers, and experts in the field. Slack is an excellent platform for real-time discussions, sharing knowledge, and getting help with FEniCS-related topics. You can join the Slack community by using [this link](https://join.slack.com/t/fenicsworkshop/shared_invite/zt-1zqrdlvmr-LziJQ59NdEOBzn7YZf59hg). 
4. **FEniCS online documentation:** The online documentation since V1.3.0 is available at [this link](https://fenicsproject.org/olddocs/dolfin/). You can go into the python folder of any version to find the syntax and meaning of different commands. The latest version is not that descriptive, and thus I recommend checking out older documentation, as most of the commands are the same.

Note: Please go through the official community [page here](https://fenicsproject.org/community/).