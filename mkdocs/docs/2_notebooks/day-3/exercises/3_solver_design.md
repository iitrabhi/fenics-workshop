# Solver design

The success of tackling complex problems often relies on the effectiveness of the solver employed. While formulating the mathematical model is crucial, paying attention to the design of an efficient and accurate solver is equally vital, as it directly impacts the quality and reliability of the results you obtain.

In this exercise you will study the impact of different solver parameters on the solution.


```python
from dolfin import *

def elasticity_problem(num_ele_along_depth=30):
    length, depth = 3, .300
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
    print("Number of degree's of freedom {}".format(U.dim()))
    problem = LinearVariationalProblem(a, L, u, bc)
    return problem, u
```

https://fenicsproject.org/pub/tutorial/html/._ftut1017.html
|Solver | Description | | | Preconditioner    |Description   |
| -- |--|--|--|--|--|
| bicgstab       |  Biconjugate gradient stabilized method| | | amg              |  Algebraic multigrid|
| cg             |  Conjugate gradient method| | | default          |  default preconditioner|
| default        |  default linear solver| | | hypre_amg        |  Hypre algebraic multigrid (BoomerAMG)|
| gmres          |  Generalized minimal residual method| | | hypre_euclid     |  Hypre parallel incomplete LU factorization|
| minres         |  Minimal residual method| | | hypre_parasails  |  Hypre parallel sparse approximate inverse|
| mumps          |  MUMPS (MUltifrontal Massively Parallel Sparse direct Solver)| | | icc              |  Incomplete Cholesky factorization|
| petsc          |  PETSc built in LU solver| | | ilu              |  Incomplete LU factorization|
| richardson     |  Richardson method| | | jacobi           |  Jacobi iteration|
| superlu        |  SuperLU| | | none             |  No preconditioner|
| tfqmr          |  Transpose-free quasi-minimal residual method| | | petsc_amg        |  PETSc algebraic multigrid|
| umfpack        |  UMFPACK (Unsymmetric MultiFrontal sparse LU factorization)| | | sor              |  Successive over-relaxation|




```python
problem, u = elasticity_problem(num_ele_along_depth=30)
solver = LinearVariationalSolver(problem)

prm = solver.parameters
prm['linear_solver'] = 'cg'
prm['preconditioner'] = 'ilu'
prm['krylov_solver']['absolute_tolerance'] = 1E-9
prm['krylov_solver']['relative_tolerance'] = 1E-9
prm['krylov_solver']['maximum_iterations'] = 1000

solver.solve()
print("The minimum displacement is: {0:6.3e} m".format(u.vector().min()))
```

    Number of degree's of freedom 18662
    The minimum displacement is: -4.712e-04 m


## The minimum displacement is: -4.71e-04 m
- For each parameter variation, record the solver parameters used and the corresponding solution time and differnce in the minimum displacement.
- Compare the results obtained with different solver parameters.
- Analyze how the solution time and accuracy are affected by varying the solver parameters.
- Based on your observations, discuss which solver parameter(s) seem to have the most significant impact on solution accuracy and computational time.
- Consider the trade-offs between accuracy and computation time when choosing different solver configurations.

1. Task 1: Considering `num_ele_along_depth=30`, change the tolerance to 1E-5. What difference do you observe in the solution?
2. Task 2: Now change the preconditioner to `hypre_euclid`. What happens to the solution?
3. Task 3: Change the preconditioner to `none` and see what happens
4. Task 4: Change the solver to `mumps` and see what happens


