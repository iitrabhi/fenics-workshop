## Plane stress v/s plane strain

The objective of this exercise is to implement a finite element analysis (FEA) program using FEniCS to simulate the behavior of a beam under bending, considering both plane stress and plane strain conditions. You will write a FEniCS code to perform the analysis for both cases and compare the results with a 3D analysis of the same beam.

<?xml version="1.0" encoding="UTF-8" standalone="no"?><!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd"><svg width="80%" height="80%" viewBox="0 0 2520 690" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xml:space="preserve" xmlns:serif="http://www.serif.com/" style="fill-rule:evenodd;clip-rule:evenodd;stroke-linejoin:round;stroke-miterlimit:2;"><rect id="Artboard1" x="0" y="0" width="2519.93" height="689.376" style="fill:#fff;"/><rect x="120.44" y="128.389" width="1514.81" height="450" style="fill:#b9b9b9;"/><rect x="84.62" y="120.54" width="35.82" height="465.699" style="fill:#3f4a47;"/><rect x="1870.64" y="119.688" width="300" height="450" style="fill:#de8383;"/><text x="813.67px" y="98.957px" style="font-family:'Montserrat-Bold', 'Montserrat';font-weight:700;font-size:93.01px;fill:#3b3b3b;">1m</text><text x="2192.32px" y="377.242px" style="font-family:'Montserrat-Bold', 'Montserrat';font-weight:700;font-size:93.01px;fill:#3b3b3b;">0<tspan x="2254.26px " y="377.242px ">.</tspan>45m</text><text x="1902.06px" y="652.461px" style="font-family:'Montserrat-Bold', 'Montserrat';font-weight:700;font-size:93.01px;fill:#3b3b3b;">0<tspan x="1964px 1989.3px " y="652.461px 652.461px ">.3</tspan>m</text><path d="M1689.89,555.343l-10.318,10.318c-0,-0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,-0c-1.756,1.756 -1.756,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,-0Z"/><path d="M1675.07,128.389l-0,443.636c-0,2.484 2.016,4.5 4.5,4.5c2.483,-0 4.5,-2.016 4.5,-4.5l-0,-443.636c-0,-2.484 -2.017,-4.5 -4.5,-4.5c-2.484,-0 -4.5,2.016 -4.5,4.5Z"/><path d="M1916,171.816l-10.318,10.318c0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.607,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.757,-1.756 1.757,-4.608 0,-6.364c-1.756,-1.756 -4.607,-1.756 -6.364,0Z"/><path d="M1901.18,119.688l0,68.81c0,2.484 2.017,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l0,-68.81c0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.483,0 -4.5,2.016 -4.5,4.5Z"/><path d="M1916,279.067l-10.318,10.318c0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.607,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.757,-1.756 1.757,-4.608 0,-6.364c-1.756,-1.756 -4.607,-1.756 -6.364,0Z"/><path d="M1901.18,226.939l0,68.81c0,2.484 2.017,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l0,-68.81c0,-2.483 -2.016,-4.5 -4.5,-4.5c-2.483,0 -4.5,2.017 -4.5,4.5Z"/><path d="M1916,396.816l-10.318,10.318c0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.607,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.757,-1.756 1.757,-4.608 0,-6.364c-1.756,-1.756 -4.607,-1.756 -6.364,0Z"/><path d="M1901.18,344.688l0,68.81c0,2.484 2.017,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l0,-68.81c0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.483,0 -4.5,2.016 -4.5,4.5Z"/><path d="M1916,508.49l-10.318,10.318c0,-0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.607,-1.756 -6.364,-0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.757,-1.756 1.757,-4.608 0,-6.364c-1.756,-1.756 -4.607,-1.756 -6.364,-0Z"/><path d="M1901.18,456.362l0,68.81c0,2.483 2.017,4.5 4.5,4.5c2.484,-0 4.5,-2.017 4.5,-4.5l0,-68.81c0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.483,-0 -4.5,2.016 -4.5,4.5Z"/><path d="M1960.32,171.816l-10.318,10.318c-0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M1945.5,119.688l-0,68.81c-0,2.484 2.016,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l-0,-68.81c-0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.484,0 -4.5,2.016 -4.5,4.5Z"/><path d="M1960.32,279.067l-10.318,10.318c-0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M1945.5,226.939l-0,68.81c-0,2.484 2.016,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l-0,-68.81c-0,-2.483 -2.016,-4.5 -4.5,-4.5c-2.484,0 -4.5,2.017 -4.5,4.5Z"/><path d="M1960.32,396.816l-10.318,10.318c-0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M1945.5,344.688l-0,68.81c-0,2.484 2.016,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l-0,-68.81c-0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.484,0 -4.5,2.016 -4.5,4.5Z"/><path d="M1960.32,508.49l-10.318,10.318c-0,-0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,-0c-1.756,1.756 -1.756,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,-0Z"/><path d="M1945.5,456.362l-0,68.81c-0,2.483 2.016,4.5 4.5,4.5c2.484,-0 4.5,-2.017 4.5,-4.5l-0,-68.81c-0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.484,-0 -4.5,2.016 -4.5,4.5Z"/><path d="M2004.64,171.816l-10.318,10.318c0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.607,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.757,-1.756 1.757,-4.608 0,-6.364c-1.756,-1.756 -4.607,-1.756 -6.364,0Z"/><path d="M1989.82,119.688l0,68.81c0,2.484 2.017,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l0,-68.81c0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.483,0 -4.5,2.016 -4.5,4.5Z"/><path d="M2004.64,279.067l-10.318,10.318c0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.607,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.757,-1.756 1.757,-4.608 0,-6.364c-1.756,-1.756 -4.607,-1.756 -6.364,0Z"/><path d="M1989.82,226.939l0,68.81c0,2.484 2.017,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l0,-68.81c0,-2.483 -2.016,-4.5 -4.5,-4.5c-2.483,0 -4.5,2.017 -4.5,4.5Z"/><path d="M2004.64,396.816l-10.318,10.318c0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.607,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.757,-1.756 1.757,-4.608 0,-6.364c-1.756,-1.756 -4.607,-1.756 -6.364,0Z"/><path d="M1989.82,344.688l0,68.81c0,2.484 2.017,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l0,-68.81c0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.483,0 -4.5,2.016 -4.5,4.5Z"/><path d="M2004.64,508.49l-10.318,10.318c0,-0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.607,-1.756 -6.364,-0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.757,-1.756 1.757,-4.608 0,-6.364c-1.756,-1.756 -4.607,-1.756 -6.364,-0Z"/><path d="M1989.82,456.362l0,68.81c0,2.483 2.017,4.5 4.5,4.5c2.484,-0 4.5,-2.017 4.5,-4.5l0,-68.81c0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.483,-0 -4.5,2.016 -4.5,4.5Z"/><path d="M2048.96,171.816l-10.318,10.318c-0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M2034.14,119.688l-0,68.81c-0,2.484 2.016,4.5 4.5,4.5c2.483,0 4.5,-2.016 4.5,-4.5l-0,-68.81c-0,-2.484 -2.017,-4.5 -4.5,-4.5c-2.484,0 -4.5,2.016 -4.5,4.5Z"/><path d="M2048.96,279.067l-10.318,10.318c-0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M2034.14,226.939l-0,68.81c-0,2.484 2.016,4.5 4.5,4.5c2.483,0 4.5,-2.016 4.5,-4.5l-0,-68.81c-0,-2.483 -2.017,-4.5 -4.5,-4.5c-2.484,0 -4.5,2.017 -4.5,4.5Z"/><path d="M2048.96,396.816l-10.318,10.318c-0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M2034.14,344.688l-0,68.81c-0,2.484 2.016,4.5 4.5,4.5c2.483,0 4.5,-2.016 4.5,-4.5l-0,-68.81c-0,-2.484 -2.017,-4.5 -4.5,-4.5c-2.484,0 -4.5,2.016 -4.5,4.5Z"/><path d="M2048.96,508.49l-10.318,10.318c-0,-0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,-0c-1.756,1.756 -1.756,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,-0Z"/><path d="M2034.14,456.362l-0,68.81c-0,2.483 2.016,4.5 4.5,4.5c2.483,-0 4.5,-2.017 4.5,-4.5l-0,-68.81c-0,-2.484 -2.017,-4.5 -4.5,-4.5c-2.484,-0 -4.5,2.016 -4.5,4.5Z"/><path d="M2093.28,171.816l-10.318,10.318c0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M2078.46,119.688l0,68.81c0,2.484 2.017,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l0,-68.81c0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.483,0 -4.5,2.016 -4.5,4.5Z"/><path d="M2093.28,279.067l-10.318,10.318c0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M2078.46,226.939l0,68.81c0,2.484 2.017,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l0,-68.81c0,-2.483 -2.016,-4.5 -4.5,-4.5c-2.483,0 -4.5,2.017 -4.5,4.5Z"/><path d="M2093.28,396.816l-10.318,10.318c0,0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M2078.46,344.688l0,68.81c0,2.484 2.017,4.5 4.5,4.5c2.484,0 4.5,-2.016 4.5,-4.5l0,-68.81c0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.483,0 -4.5,2.016 -4.5,4.5Z"/><path d="M2093.28,508.49l-10.318,10.318c0,-0 -10.318,-10.318 -10.318,-10.318c-1.756,-1.756 -4.608,-1.756 -6.364,-0c-1.756,1.756 -1.756,4.608 0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 0,-6.364c-1.756,-1.756 -4.608,-1.756 -6.364,-0Z"/><path d="M2078.46,456.362l0,68.81c0,2.483 2.017,4.5 4.5,4.5c2.484,-0 4.5,-2.017 4.5,-4.5l0,-68.81c0,-2.484 -2.016,-4.5 -4.5,-4.5c-2.483,-0 -4.5,2.016 -4.5,4.5Z"/><path d="M2137.6,171.816l-10.318,10.318c-0,0 -10.318,-10.318 -10.318,-10.318c-1.757,-1.756 -4.608,-1.756 -6.364,0c-1.757,1.756 -1.757,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.757,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M2122.78,119.688l-0,68.81c-0,2.484 2.016,4.5 4.5,4.5c2.483,0 4.5,-2.016 4.5,-4.5l-0,-68.81c-0,-2.484 -2.017,-4.5 -4.5,-4.5c-2.484,0 -4.5,2.016 -4.5,4.5Z"/><path d="M2137.6,279.067l-10.318,10.318c-0,0 -10.318,-10.318 -10.318,-10.318c-1.757,-1.756 -4.608,-1.756 -6.364,0c-1.757,1.756 -1.757,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.757,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M2122.78,226.939l-0,68.81c-0,2.484 2.016,4.5 4.5,4.5c2.483,0 4.5,-2.016 4.5,-4.5l-0,-68.81c-0,-2.483 -2.017,-4.5 -4.5,-4.5c-2.484,0 -4.5,2.017 -4.5,4.5Z"/><path d="M2137.6,396.816l-10.318,10.318c-0,0 -10.318,-10.318 -10.318,-10.318c-1.757,-1.756 -4.608,-1.756 -6.364,0c-1.757,1.756 -1.757,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.757,-1.756 -4.608,-1.756 -6.364,0Z"/><path d="M2122.78,344.688l-0,68.81c-0,2.484 2.016,4.5 4.5,4.5c2.483,0 4.5,-2.016 4.5,-4.5l-0,-68.81c-0,-2.484 -2.017,-4.5 -4.5,-4.5c-2.484,0 -4.5,2.016 -4.5,4.5Z"/><path d="M2137.6,508.49l-10.318,10.318c-0,-0 -10.318,-10.318 -10.318,-10.318c-1.757,-1.756 -4.608,-1.756 -6.364,-0c-1.757,1.756 -1.757,4.608 -0,6.364l16.682,16.682l16.682,-16.682c1.756,-1.756 1.756,-4.608 -0,-6.364c-1.757,-1.756 -4.608,-1.756 -6.364,-0Z"/><path d="M2122.78,456.362l-0,68.81c-0,2.483 2.016,4.5 4.5,4.5c2.483,-0 4.5,-2.017 4.5,-4.5l-0,-68.81c-0,-2.484 -2.017,-4.5 -4.5,-4.5c-2.484,-0 -4.5,2.016 -4.5,4.5Z"/></svg>

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
