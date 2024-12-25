The accompanying Jupyter notebook can be obtained here [1_bi_metallic_quasi_static](../../../../../src/day-5/exercises/1_bi_metallic_quasi_static.ipynb)

**Exercise: Bimetallic Quasi-Static Thermal Simulation using Steel and Brass Unit Square using FEniCS**

In this exercise, you will perform a bimetallic quasi-static thermal simulation using the Finite Element Method (FEM) with FEniCS. The simulation will involve a unit square made of two different materials: steel and brass. The goal is to analyze the steady-state temperature distribution within the square and observe the thermal behavior at the interface of the two materials.     |



**Materials:**

The following table provides the material properties for steel and brass:

| Material | Thermal Conductivity (k) [W/mK] | Specific Heat Capacity (Cp) [J/kgK] | Density (ρ) [kg/m³] |
|----------|-----------------------------|-----------------------------------|--------------------|
| Steel    | 50.2                        | 486.0                             | 7850               |
| Brass    | 109.0                       | 377.0                             | 8520               |


**Task:**

1. **Geometry and Mesh Generation:**
   - Create a unit square with side length L = 1.0 units (you can choose any unit system).
   - Divide the square into a suitable number of elements to create a mesh. You can start with a relatively coarse mesh and later refine it to observe its effect on the simulation results.



2. **Setting up the Problem:**
   - Define the governing equation for steady-state heat conduction in 2D. The equation involves the Laplace operator and accounts for the material properties.
   - Implement appropriate boundary conditions for the simulation. Specify the temperature boundary condition for all four sides of the unit square. You can use Dirichlet boundary conditions, where you set the temperature values at the boundaries.



3. **Material Properties:**
   - Use the provided material properties for steel and brass.



4. **Thermal Simulation:**
   - Assemble the finite element problem using FEniCS and solve the system of equations to obtain the temperature distribution within the unit square.
   - Perform the simulation for different time steps to observe the quasi-static behavior of the system. You can use a small time step for better accuracy.



5. **Post-Processing and Analysis:**
   - Visualize the temperature distribution using contour plots or color maps.
   - Analyze the temperature distribution at the interface between steel and brass. Observe if there are any significant temperature gradients or discontinuities at this interface.



6. **Mesh Refinement:**
   - Re-run the simulation with a refined mesh (i.e., more elements) and compare the results with the coarser mesh. Discuss the differences and the effect of mesh refinement on the accuracy of the simulation.



7. **Discussion and Conclusion:**
   - Write a summary of your findings, including observations about the temperature distribution, behavior at the interface of the two materials, and the influence of mesh refinement on the simulation results.

