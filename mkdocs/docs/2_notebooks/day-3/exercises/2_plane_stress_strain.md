The accompanying Jupyter notebook can be obtained here [2_plane_stress_strain](../../../src/day-3/exercises/2_plane_stress_strain.ipynb)



## Plane stress v/s plane strain

The objective of this exercise is to implement a finite element analysis (FEA) program using FEniCS to simulate the behavior of a beam under bending, considering both plane stress and plane strain conditions. You will write a FEniCS code to perform the analysis for both cases and compare the results with a 3D analysis of the same beam.

![png](2_plane_stress_strain_files/2_plane_stress_strain_image.png)

Materials:
- Steel material properties:
  - Young's Modulus (E): 200 GPa
  - Poisson's Ratio (ν): 0.3
  - Density (ρ): 7850 kg/m^3

Bar Dimensions:
- Length (Lx): 1 meter
- Ly         : 0.3 meter 
- Lz         : 0.45 meter

Boundary Conditions:
- One end of the beam is fixed (fixed boundary condition).
- The other end is subjected to a uniaxial downward load:
  - Load (F): 10,000 N

Mesh:
- For 3D analysis use a simple 3D mesh for the beam, consisting of tetrahedral elements.
- For 2D analysis use a simple 2D mesh for the beam, consisting of triangular elements.

---
Steps:
- Perform 3D analysis
    1. Import FEniCS and Necessary Libraries
    2. Define the 3D Geometry and Mesh
    3. Define the Material Properties
    4. Define the Boundary Conditions
    5. Define the Finite Element Function Space
    6. Formulate the Governing Equations
    7. Solve the System
- Perform 2D plane stress analysis
    1. Define the 2D Geometry and Mesh
    2. Define the Material Properties
    3. Define the Boundary Conditions
    4. Define the Finite Element Function Space
    5. Formulate the Governing Equations
    6. Solve the System
- Perform 2D plane strain analysis
    1. Define the 2D Geometry and Mesh
    2. Define the Material Properties
    3. Define the Boundary Conditions
    4. Define the Finite Element Function Space
    5. Formulate the Governing Equations
    6. Solve the System
- Compare the results of plane stress, plane strain and 3D analysis with analytical solution.
- See what happens to the solution when you increase the mesh density


```python
sim_types = ["3d","plane_stress", "plane_strain"]
sim_type = sim_types[0]

E0, nu = 2e11, 0.3
mu = E0 / (2 * (1 + nu))
lmbda = E0 * nu / ((1 + nu) * (1 - 2 * nu))

if sim_type == "plane_stress":
    lmbda = 2 * mu * lmbda / (lmbda + 2 * mu)
```


```python

```
