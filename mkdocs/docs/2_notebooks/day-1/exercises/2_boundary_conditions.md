# Boundary conditions

In FEniCS, the "CompiledSubDomain" class is a useful tool that allows users to define complex subdomains within a given computational domain for finite element simulations. Subdomains are portions of the computational domain where different physical or material properties are applied, or specific boundary conditions are imposed.

The primary advantage of using the "CompiledSubDomain" class is that it allows you to define subdomains using mathematical expressions or conditions, which are then compiled into efficient C++ code. This compiled code is utilized during the simulation, providing a significant performance boost compared to interpreting the subdomain expressions directly in Python.

Please visit the official documentation link provided to learn how to modify the bounday conditions using "CompiledSubDomain". After familiarizing yourself with the process, return here to implement the changes and visualize the updated results.

https://fenicsproject.org/pub/tutorial/sphinx1/._ftut1005.html#using-c-code-snippets-to-define-subdomains


```python
from dolfin import *
```

**Change this portion of the code to:**
1. Make a 2D unit square mesh.


```python
mesh = IntervalMesh(30, 0, 1)
```


```python
U = FunctionSpace(mesh, "CG", 1)
```

**Change this portion of the code to:**
1. Mark only the left edge as fixed.
2. Mark the left and top edge as fixed.


```python
u_D = Constant(0.0)
boundary = CompiledSubDomain("on_boundary")
bc = DirichletBC(U, u_D, boundary)
```


```python
u, v = TrialFunction(U), TestFunction(U)

a = inner(grad(u), grad(v)) * dx
f_expr = Expression("pi*pi*sin(pi*x[0])", pi=np.pi, degree=2)
L = f_expr * v * dx

u_sol = Function(U, name = "field")
solve(a == L, u_sol, bc)

with XDMFFile("output/result.xdmf") as outfile:
    outfile.write(u_sol)
```
