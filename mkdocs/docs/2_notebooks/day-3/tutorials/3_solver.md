The accompanying Jupyter notebook can be obtained here [3_solver](../../../src/day-3/tutorials/3_solver.ipynb)



## Solver Design

In addition to the standard `solve` method, FEniCS offers a powerful tool called the `LinearVariationalSolver`, which grants users the ability to finely adjust and control various parameters of the `solver`. This enhanced control allows for precise customization and optimization of the solver's behavior, leading to improved accuracy and efficiency in solving partial differential equations.

In this tutorial you will learn how to define a Linear variational problem and modify the solver parameters.



```python
from dolfin import *
length, depth = 3, .300
num_ele_along_depth = 300
ele_size = depth/num_ele_along_depth
mesh = RectangleMesh(Point(0, 0), Point(length, depth),
                     int(length/ele_size), int(depth/ele_size))
U = VectorFunctionSpace(mesh, 'CG', 1)
dim = mesh.topology().dim()
clamped_boundary = CompiledSubDomain("near(x[0],0)")
bc = DirichletBC(U, Constant((0,)*dim), clamped_boundary)
E, nu = 2e11, 0.3
rho, g = 7800, 9.81
lmbda = (E * nu) / ((1 + nu) * (1 - 2 * nu))
mu = E / (2 * (1 + nu))


def epsilon(u):
    return 0.5*(grad(u) + grad(u).T)


def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(dim) + 2*mu*epsilon(u)


# Define variational problem
u, v = TrialFunction(U), TestFunction(U)
f = Constant((0, -rho*g))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx

u = Function(U)
```


```python
print("Number of degree's of freedom {}".format(U.dim()))
```

    Number of degree's of freedom 1806602



```python
# Compute solution
solve(a == L, u, bc)
u.vector().min()
```


    ---------------------------------------------------------------------------

    RuntimeError                              Traceback (most recent call last)

    <ipython-input-118-fed64cb3fd67> in <module>
          1 # Compute solution
          2 u = Function(U)
    ----> 3 solve(a == L, u, bc)
          4 u.vector().min()


    /usr/local/lib/python3.6/dist-packages/dolfin/fem/solving.py in solve(*args, **kwargs)
        218     # tolerance)
        219     elif isinstance(args[0], ufl.classes.Equation):
    --> 220         _solve_varproblem(*args, **kwargs)
        221 
        222     # Default case, just call the wrapped C++ solve function


    /usr/local/lib/python3.6/dist-packages/dolfin/fem/solving.py in _solve_varproblem(*args, **kwargs)
        245         solver = LinearVariationalSolver(problem)
        246         solver.parameters.update(solver_parameters)
    --> 247         solver.solve()
        248 
        249     # Solve nonlinear variational problem


    RuntimeError: 
    
    *** -------------------------------------------------------------------------
    *** DOLFIN encountered an error. If you are not able to resolve this issue
    *** using the information listed below, you can ask for help at
    ***
    ***     fenics-support@googlegroups.com
    ***
    *** Remember to include the error message listed below and, if possible,
    *** include a *minimal* running example to reproduce the error.
    ***
    *** -------------------------------------------------------------------------
    *** Error:   Unable to successfully call PETSc function 'KSPSolve'.
    *** Reason:  PETSc error code is: 76 (Error in external library).
    *** Where:   This error was encountered inside /tmp/dolfin/dolfin/la/PETScKrylovSolver.cpp.
    *** Process: 0
    *** 
    *** DOLFIN version: 2019.1.0
    *** Git changeset:  74d7efe1e84d65e9433fd96c50f1d278fa3e3f3f
    *** -------------------------------------------------------------------------



 FEniCS' standard solve method relies on a direct solver, which proves inadequate for computing solutions in systems with degrees of freedom exceeding approximately one million. To address this limitation, iterative solvers and preconditioners become necessary alternatives to efficiently handle large-scale problems. By employing these techniques, FEniCS enables the successful computation of solutions in scenarios where the direct solver falls short, making it a valuable tool for tackling complex simulations and high-dimensional models.


```python
problem = LinearVariationalProblem(a, L, u, bc)
solver = LinearVariationalSolver(problem)

prm = solver.parameters
prm['linear_solver'] = 'cg'
prm['preconditioner'] = 'hypre_euclid'
prm['krylov_solver']['absolute_tolerance'] = 1E-5
prm['krylov_solver']['relative_tolerance'] = 1E-5
prm['krylov_solver']['maximum_iterations'] = 1000

solver.solve()
```


```python
print("The minimum displacement is: {0:6.3e} m".format(u.vector().min()))
```

    The minimum displacement is: -4.733e-04 m



```python

```
