## Spacially varying properties

With FEniCS, researchers and engineers can model and simulate problems with spatially varying material properties efficiently and accurately. We can use the Expression class provided by FEniCS to represent the material property as a function of the co-ordinate of the mesh.

Please visit the official documentation link provided to learn how to modify expressions. After familiarizing yourself with the process, return here to implement the changes and visualize the updated results.

<https://fenicsproject.org/pub/tutorial/sphinx1/._ftut1003.html#index-28>


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
1. Change the material property to $E(x) = 10 + x$


```python
E = Constant(10)
```


```python
u_D = Constant(0.0)
boundary = CompiledSubDomain("on_boundary")
bc = DirichletBC(U, u_D, boundary)
u, v = TrialFunction(U), TestFunction(U)

a = E * inner(grad(u), grad(v)) * dx
f_expr = Expression("pi*pi*sin(pi*x[0])", pi=np.pi, degree=2)
L = f_expr * v * dx

u_sol = Function(U, name = "field")
solve(a == L, u_sol, bc)

with XDMFFile("output/result.xdmf") as outfile:
    outfile.write(u_sol)
```

    Calling FFC just-in-time (JIT) compiler, this may take some time.



```python

```
