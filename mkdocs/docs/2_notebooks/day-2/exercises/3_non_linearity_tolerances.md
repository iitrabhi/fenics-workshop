## Tolerances and convergence

In numerical methods, solving non-linear problems poses unique challenges due to their intrinsic complexity. Solver design become crucial aspect to consider while tackling these problems, as they govern the accuracy and reliability of the obtained solutions. 

Understanding tolerance levels, which dictate the acceptable deviation between successive iterations, and convergence criteria, which indicate when a solution has reached a satisfactory result is of prime importance when handling non-linear problems.

In this exercise, your task is to explore the impact of modifying three parameters: absolute tolerance, relative tolerance, and maximum iterations, on the computational cost and solution accuracy. By adjusting these parameters, you will gain insights into how they influence the numerical solution of the problem at hand.


```
from dolfin import *
import numpy

mesh = IntervalMesh(40, 0, 1)
V = FunctionSpace(mesh, 'Lagrange', 1)

left_boundary = CompiledSubDomain("on_boundary && near(x[0],0)")
right_boundary = CompiledSubDomain("on_boundary && near(x[0],1)")

bc_0 = DirichletBC(V, Constant(0.0), left_boundary)
bc_1 = DirichletBC(V, Constant(1.0), right_boundary)
bcs = [bc_0, bc_1]

m = 5


def q(u):
    return (1+u)**m

# Define variational problem
v = TestFunction(V)
du = TrialFunction(V)
u = Function(V)  # most recently computed solution
F = inner(q(u)*nabla_grad(u), nabla_grad(v))*dx
```

Change the tangent modifer to a value between 0.5 and 1 and see how it affects the solution.


```
tangent_modifier = 1
J = derivative(tangent_modifier*F, u, du)
```


```
# Compute solution
problem = NonlinearVariationalProblem(F, u, bcs, J)
solver = NonlinearVariationalSolver(problem)
```

Change the tolerances and maximum_iterations to see what happens to the solution in terms of computational cost and solution accuracy.


```
prm = solver.parameters
prm['newton_solver']['absolute_tolerance'] = 1E-5
prm['newton_solver']['relative_tolerance'] = 1E-5
prm['newton_solver']['maximum_iterations'] = 25

iterations, converged = solver.solve()
print("Number of iterations: {}".format(iterations))
print("The solver converged." if converged else "The solver did not converge.")
```

    Number of iterations: 6
    The solver converged.



```
# Find max error
u_exact = Expression(
    'pow((pow(2, m+1)-1)*x[0] + 1, 1.0/(m+1)) - 1', m=m, degree=1)
u_e = interpolate(u_exact, V)
diff = numpy.abs(u_e.vector()[:] - u.vector()[:]).max()
print('Max error:{0:5.3e}'.format(diff))
```

    Max error:1.559e-06



```

```


```

```
