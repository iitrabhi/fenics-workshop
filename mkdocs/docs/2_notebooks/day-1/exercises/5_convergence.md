The accompanying Jupyter notebook can be obtained here [5_convergence](../../../../../src/day-1/exercises/5_convergence.ipynb)

## Convergence

Convergence analysis in the FEM is a critical process that evaluates the accuracy and reliability of numerical solutions to PDEs. It involves discretizing the continuous problem into a mesh of elements and studying how the numerical solution converges to the exact solution as the mesh is refined. 

By analyzing various error norms and convergence rates, researchers can determine the reliability of the numerical scheme, select appropriate mesh resolutions, and validate results against analytical solutions. 

Please visit the official documentation link provided to learn how to modify the mesh. After familiarizing yourself with the process, return here to implement the changes and visualize the updated results.

<https://fenicsproject.org/olddocs/dolfin/latest/python/demos/built-in-meshes/demo_built-in-meshes.py.html>


```python
from dolfin import *
```

**Change this portion of the code to make a 2D unit square mesh and perform the analysis with mesh size:**
1. $2\times 2$
2. $5\times 5$
3. $10\times 10$
4. $20\times 20$
5. $50\times 50$


```python
mesh = IntervalMesh(30, 0, 1)
```


```python
U = FunctionSpace(mesh, "CG", 1)

u_D = Constant(0.0)
boundary = CompiledSubDomain("on_boundary")
bc = DirichletBC(U, u_D, boundary)
u, v = TrialFunction(U), TestFunction(U)

a = inner(grad(u), grad(v)) * dx
f_expr = Expression("pi*pi*sin(pi*x[0])", pi=np.pi, degree=2)
L = f_expr * v * dx

u_sol = Function(U, name = "field")
solve(a == L, u_sol, bc)

with XDMFFile("output/result.xdmf") as outfile:
    outfile.write(u_sol)
```
