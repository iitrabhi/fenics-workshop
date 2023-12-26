**Thermo-Mechanical Transient Simulation using FEniCS**

**Objective:**
The objective of this exercise is to perform a thermo-mechanical transient simulation using the FEniCS finite element library. In this exercise, we will explore the coupling between thermal and mechanical phenomena in a time-dependent setting. The problem involves a 2D domain, and we will consider a linear elastic material model.

**Problem Description:**
Consider a 2D square domain Ω with dimensions LxL, where L = 1.0 m. The domain is initially at room temperature T_initial = 25°C. At t = 0, the left boundary (x = 0) is subjected to a temperature of T_boundary = 100°C, while the right boundary (x = L) is maintained at T_boundary = 50°C. The top (y = L) and bottom (y = 0) boundaries are assumed to be insulated, i.e., there is no heat flux through these boundaries.

The material properties are given in the following table:

| Property         | Symbol   | Value  | 
|------------------|----------|---------------------------|
| Young's modulus  | E        | 210e9 Pa                  | 
| Poisson's ratio  | ν        | 0.3                       | 
| Thermal expansion| α        | 1.2e-5 1/°C              | 
| Conductivity     | k        | 50 W/(m°C)                | 
| Density          | ρ        | 7800 kg/m³                | 
| Specific heat    | c        | 480 J/(kg°C)              |


**Tasks:**
1. Set up the 2D square domain Ω and define appropriate mesh resolution.
2. Implement a function to compute the temperature distribution within the domain over time using the heat equation.
3. Implement a function to compute the displacement field within the domain over time using the linear elasticity equations.
4. Perform a transient simulation, coupling the heat equation and linear elasticity equations, to obtain the temperature and displacement fields at different time steps.
5. Visualize the temperature and displacement fields at various time steps using appropriate visualization tools.
6. Analyze the results and observe the temperature and displacement evolution within the domain over time.

**Boundary Conditions:**
- Initial condition: T(x, y, 0) = T_initial for all points (x, y) within Ω.
- Left boundary condition: T(x = 0, y, t) = T_boundary for all y and t > 0.
- Right boundary condition: T(x = L, y, t) = T_boundary for all y and t > 0.
- Top and bottom boundary conditions: ∂T/∂y = 0 (insulated boundaries) for all x, t > 0.
- Displacement boundary conditions: u(x = 0, y, t) = 0 for all y and t > 0 (fixed boundary).

**Note:**
For the simulation, you can use an appropriate time-stepping scheme (e.g., Forward Euler, Backward Euler, or Crank-Nicolson) and consider suitable time intervals to observe the transient behavior effectively. You may also choose the finite element degree and other numerical parameters based on your understanding and computational resources.

Remember to interpret and analyze the results in light of the thermo-mechanical coupling and transient behavior exhibited by the system.
