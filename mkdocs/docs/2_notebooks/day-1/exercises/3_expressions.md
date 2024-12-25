The accompanying Jupyter notebook can be obtained here [3_expressions](../../../src/day-1/exercises/3_expressions.ipynb)



## Expressions

In FEniCS, an Expression is a flexible and convenient way to define mathematical functions or expressions within the domain of interest. Expressions are often used to specify boundary conditions, source terms, initial conditions, or any other function needed in the formulation of the PDE problem.

The beauty of using Expressions lies in their simplicity and directness. Users can define an Expression using a concise mathematical expression, incorporating spatial coordinates and/or time variables. Additionally, FEniCS supports the use of elementary mathematical functions, mathematical constants, and custom-defined functions within Expressions.

Please visit the official documentation link provided to learn how to modify expressions. After familiarizing yourself with the process, return here to implement the changes and visualize the updated results.

<https://hplgit.github.io/fenics-tutorial/pub/sphinx1/._ftut1003.html#index-28>


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

u_D = Constant(0.0)
boundary = CompiledSubDomain("on_boundary")
bc = DirichletBC(U, u_D, boundary)
u, v = TrialFunction(U), TestFunction(U)

a = inner(grad(u), grad(v)) * dx

```

**Change this portion of the code to:**
1. Change the forcing function to $f(x) = x$
2. Change the forcing function to $f(x,y) = sin(\pi x)sin(\pi y)$


```python
f_expr = Expression("pi*pi*sin(pi*x[0])", pi=np.pi, degree=2)
```


```python
L = f_expr * v * dx

u_sol = Function(U, name = "field")
solve(a == L, u_sol, bc)

with XDMFFile("output/result.xdmf") as outfile:
    outfile.write(u_sol)
```
