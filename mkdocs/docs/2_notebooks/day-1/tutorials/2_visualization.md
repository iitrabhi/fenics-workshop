The accompanying Jupyter notebook can be obtained here [2_visualization](../../../../../src/day-1/tutorials/2_visualization.ipynb)

## Visualizing simulation results

Once the simulation is performed in FEniCS and the results are obtained, this tutorial will show you how to save the output in the XDMF file format. XDMF is ideal for storing large-scale scientific data, especially for finite element simulations, thanks to its flexibility and efficiency.


```python
from dolfin import *
num_elements = 30
mesh = IntervalMesh(num_elements, 0, 1)
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
```

## Option 1: Write data efficiently using `with` syntax.

FEniCS provides us with an XDMFFile class for efficiently handling XDMF files and storing simulation results.


```python
with XDMFFile("output/result.xdmf") as outfile:
    outfile.write(u_sol)
```

Let's break down the command step by step:

1. `XDMFFile`: This is a class in FEniCS used to write simulation data to an XDMF file. XDMF is an XML-based file format commonly used to store scientific data, especially for finite element simulations.

2. `"output/result.xdmf"`: This is the file path and name where the XDMF file will be created or updated. In this case, the file will be named "result.xdmf" and will be located in the "output" directory (relative to the current working directory).

3. `with`: This keyword is used to define a context manager in Python. It ensures that resources associated with the context (in this case, the XDMFFile object) are properly managed and released when the block of code inside the `with` statement is executed.

4. `as outfile`: This assigns the XDMFFile object to the variable `outfile`, which can be used to interact with the file and write data.

5. `outfile.write(u_sol)`: This line of code writes the data `u_sol` to the XDMF file. `u_sol` is the variable containing the solution. Note that the name that is passed on to the variable during its creation is the name that will appear in the visualization tool paraview.

Overall, this command creates an XDMFFile object, opens the file "output/result.xdmf" for writing, writes the data contained in the `u_sol` variable to the file, and then automatically closes the file after the code block within the `with` statement is executed. This allows for efficient and controlled writing of simulation results to an XDMF file, which can later be visualized using tools like Paraview to analyze and interpret the simulation data.

## Option 2: Use this to write data when dealing with timeseries

In engineering, numerous problems demand the consideration of multiple variables and time-steps, which subsequently leads to the need for writing these variables to a file in a time series format. In such cases, FEniCS offers a convenient solution through the utilization of the XDMFFile class, facilitating efficient storage of simulation results.

**We can access the `parameters` of the `XDMFFile` object to control its behaviour:**
- **functions_share_mesh**: Default is false, it controls whether all functions on a single time step share the same mesh. If true the files created will be smaller and also behave better in Paraview.
- **rewrite_function_mesh**: Default settings is true, i.e, it controls whether the mesh will be rewritten every time step. If the mesh does not change this can be turned off to create smaller files. 
- **flush_output**: Default is false, it controls the ability of Paraview to render during execution. If you are doing a time dependent analysis, setting it to true will allow you to visualize results during run.


```python
outfile = XDMFFile("output/result.xdmf")
outfile.parameters["functions_share_mesh"] = True
outfile.parameters["rewrite_function_mesh"] = False
outfile.parameters["flush_output"] = True

time_step = 1
outfile.write(u_sol, time_step)
outfile.close()
```

At the end of the execution it is recommended to close the `XDMFFile` object.
