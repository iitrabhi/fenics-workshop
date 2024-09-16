# Poissons Equation

Welcome to this FEniCS tutorial, where we will explore how to verify the accuracy of a Poisson's equation solver using the "manufactured solution" technique. FEniCS is a powerful open-source finite element library for solving partial differential equations (PDEs), widely used for scientific computing and simulation.

The "manufactured solution" approach is a valuable method to validate the correctness of finite element implementations. In this technique, we first construct an exact solution to the PDE, often a smooth and analytically known function, that satisfies the given equation. Next, we compute the corresponding right-hand side of the PDE using the exact solution. By feeding the manufactured solution and the derived right-hand side into our FEniCS solver, we can compare the numerical solution with the exact solution, thus quantifying the solver's accuracy.

In this tutorial, we will focus on solving the one-dimensional Poisson's equation:

$$\begin{split}- \nabla^{2} u &= f \quad {\rm in} \ \Omega, \\
             u &= 0 \quad {\rm on} \ \Gamma_{D}, \\
             \nabla u \cdot n &= g \quad {\rm on} \ \Gamma_{N}. \\\end{split}
$$

subject to homogeneous Dirichlet boundary conditions, where u(x) is the unknown function, and f(x) is the right-hand side. We will construct a simple manufactured solution, u_exact(x), and calculate the corresponding f(x) that satisfies the equation.

Throughout the tutorial, we will cover the following steps:

- Importing the necessary modules.
- Defining the manufactured solution and its corresponding right-hand side.
- Creating the one-dimensional mesh using FEniCS.
- Defining the appropriate FunctionSpace for the problem.
- Imposing the homogeneous Dirichlet boundary conditions.
- Formulating the variational problem using FEniCS's TrialFunction and TestFunction.
- Solving the Poisson's equation using FEniCS's solve function.
- Comparing the numerical solution with the exact solution to quantify the solver's accuracy.

By the end of this tutorial, you will have a better understanding of the manufactured solution technique, its importance in validating finite element solvers, and how to implement it using FEniCS on an interval mesh. So, let's get started with our journey into the world - of FEniCS and manufactured solutions!


## Step 1: Import the necessary modules



```python
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```

## Step 2: Define the mesh


```python
# Create the mesh
num_elements = 3
# num_elements equally spaced intervals in [0, 1]
mesh = IntervalMesh(num_elements, 0, 1)
```


```python
plot(mesh)
```




    [<matplotlib.lines.Line2D at 0x7f9f10e6d4e0>]




    
![png](1_linear_poisson_files/1_linear_poisson_6_1.png)
    


## Step 3: Define the function space

The line of code `U = FunctionSpace(mesh, "CG", 1)` in FEniCS creates a function space `U` based on linear continuous elements (`CG`) defined on the mesh with degree `1` polynomial approximation. 


```python
U = FunctionSpace(mesh, "CG", 1)
```

- In above code **'CG'** stands for **Continuous Galerkin** or **Lagrange elements** continuous across element boundaries.
- **'1'** indicates that basic functions are linear i.e. the solution is approximated by piecewise linear function. If it is change to **2** the basis functions would be quadratic.
## Step 4: Define boundary condition
Please read [this article](https://www.resolvedanalytics.com/cfd/what-is-a-neumann-boundary-condition) beforehand if you have no idea about Dirichlet or Neumann boundary condition.

We create a `DirichletBC` object (`bc`) that associates the function space `U` with the boundary condition `u_D` and the subdomain `boundary`. This means that the solution function `u_sol` will have the value `0.0` on the boundary of the domain during the solution of the PDE
.


```python
u_D = Constant(0.0)
boundary = CompiledSubDomain("on_boundary")
bc = DirichletBC(U, u_D, boundary)
```

- **Constant(0.0)** just defines a constant value 
- **CompiledSubDomain** is a FEniCS class that allows the definition of a subdomain or part of the domain. It compiles the condition making it more efficient for large-scale problems or complex geometries.
- **on_boundary** is built in condition in FEniCS that returns 'True' for all points on the boundary of the domain and 'False' on inside the domain.

So, all the code does is create a subdomain names **boundary** that corresponds to all points on the boundary of the mesh. 
- In 1D boundary would consist of two end points of the interval. For $[a,b]$, the boundary will be at **a** and **b**
- In 2D however, boundaries will be edges i.e lines and it goes on like that for 3D where planes will form the boundaries.


And finally, the last line **bc = DirichletBC(U, u_D, boundary)** will apply the value of **u_D** to the **boundary** that has been defined.


## Step 5: Define weak form

Please refer [here](https://fenicsproject.org/pub/tutorial/sphinx1/._ftut1003.html) to know how the weak form was derived.

$$\begin{split}a(u, v) &= L(v) \quad \forall \ v \in V, \\
a(u, v) &= \int_{\Omega} \nabla u \cdot \nabla v \, {\rm d} x, \\
L(v)    &= \int_{\Omega} f v \, {\rm d} x.\end{split}$$


```python
u = TrialFunction(U)
v = TestFunction(U)
a = inner(grad(u), grad(v)) * dx
```

## Step 5.1: Define the manufactured rhs

For this tutorial, let's choose a simple manufactured solution. We will solve the Poisson's equation in 1D:

$$-\nabla^2 u(x) = f(x), 0 < x < 1,$$

where u(x) is the unknown function, and f(x) is the right-hand side. We will choose an analytical solution u_exact(x) that satisfies the above equation.

For this example, let's take:

$$u_{exact}(x) = sin(\pi x)$$

and calculate the corresponding f(x):

$$f(x) = -\nabla^2 u_{exact}(x) = \pi^2  sin(\pi x)$$



```python
f_expr = Expression("pi*pi*sin(pi*x[0])", pi=np.pi, degree=2)
```

### Visualize Expression 
In the given code snippet:

1. `f_val = project(f_expr, V)`: We use the `project` function to interpolate the expression `f_expr` onto the function space `V`. This creates a new function `f_val` that represents the projection of the expression `f_expr` onto the space `V`. The `project` function is useful when we want to create functions from `Expression` and visualize them.



```python
V = FunctionSpace(mesh, 'CG', 1)
f_val = project(f_expr, V)
plt.plot([f_val(x) for x in np.linspace(0,1,100)])
```




    [<matplotlib.lines.Line2D at 0x7f9ee473ad68>]




    
![png](1_linear_poisson_files/1_linear_poisson_17_1.png)
    



```python
L = f_expr * v * dx
```

## Compute the solution



```python
u_sol = Function(U)
solve(a == L, u_sol, bc)
```


```python
plot(u_sol)
```




    [<matplotlib.lines.Line2D at 0x7f9ee3061588>]




    
![png](1_linear_poisson_files/1_linear_poisson_21_1.png)
    


You can pass the co-ordinates of any point inside the mesh to get the value of any FEniCS function at that point


```python
u_sol(0.5)
```




    0.8648790643807451



## Post processing


```python
def u_exact(x):
    return np.sin(np.pi * x)
```


```python
points = np.linspace(0, 1, 100)

# Evaluate the exact solution at the mesh points
u_exact_values = np.array([u_exact(x) for x in points])

# Evaluate the numerical solution at the mesh points
u_numerical_values = np.array([u_sol(x) for x in points])

# Compute the error
error = u_exact_values - u_numerical_values

print("L2 error:", np.linalg.norm(error))
```

    L2 error: 0.6951127497810745



```python
plt.figure()
plt.plot(points, u_exact_values, "--", label='Exact solution', linewidth=3)
plt.plot(points, u_numerical_values, label='Numerical solution')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend()
plt.title('Manufactured solution of Poisson\'s equation')
plt.grid()
plt.show()
```


    
![png](1_linear_poisson_files/1_linear_poisson_27_0.png)
    

```
!pip install plotly
```

```python
import plotly.express as px

# Convert mesh coordinates and solutions to DataFrame for Plotly Express
import pandas as pd
data = pd.DataFrame({'x': points,
                     'Exact solution': [u_exact(x) for x in points],
                     'Numerical solution': [u_sol(x) for x in points]})

# Create a Plotly Express figure with colors specified
fig = px.line(data, x='x', y=['Exact solution', 'Numerical solution'], title='Manufactured solution of Poisson\'s equation',
              color_discrete_map={'Exact solution': 'black', 'Numerical solution': 'blue'})

fig.update_traces(mode="lines", hovertemplate=None)

# Set the figure size to achieve a 1:1 aspect ratio
fig.update_layout(
    hovermode='x',
    autosize=False,
    width=800,  # You can adjust this value to get the desired aspect ratio
    height=500
)

# Show the Plotly Express figure
fig.show()
```

If you want to save the image

```
!pip install kaleido
```

```python
fig.write_image("img.png") 
```


<div>                            <div id="9795c319-ae78-4999-8578-d00e1dfa8c61" class="plotly-graph-div" style="height:500px; width:800px;"></div>            <script type="text/javascript">                require(["plotly"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("9795c319-ae78-4999-8578-d00e1dfa8c61")) {                    Plotly.newPlot(                        "9795c319-ae78-4999-8578-d00e1dfa8c61",                        [{"legendgroup":"Exact solution","line":{"color":"black","dash":"solid"},"marker":{"symbol":"circle"},"mode":"lines","name":"Exact solution","orientation":"v","showlegend":true,"x":[0.0,0.010101010101010102,0.020202020202020204,0.030303030303030304,0.04040404040404041,0.05050505050505051,0.06060606060606061,0.07070707070707072,0.08080808080808081,0.09090909090909091,0.10101010101010102,0.11111111111111112,0.12121212121212122,0.13131313131313133,0.14141414141414144,0.15151515151515152,0.16161616161616163,0.17171717171717174,0.18181818181818182,0.19191919191919193,0.20202020202020204,0.21212121212121213,0.22222222222222224,0.23232323232323235,0.24242424242424243,0.25252525252525254,0.26262626262626265,0.27272727272727276,0.2828282828282829,0.29292929292929293,0.30303030303030304,0.31313131313131315,0.32323232323232326,0.33333333333333337,0.3434343434343435,0.3535353535353536,0.36363636363636365,0.37373737373737376,0.38383838383838387,0.393939393939394,0.4040404040404041,0.4141414141414142,0.42424242424242425,0.43434343434343436,0.4444444444444445,0.4545454545454546,0.4646464646464647,0.4747474747474748,0.48484848484848486,0.494949494949495,0.5050505050505051,0.5151515151515152,0.5252525252525253,0.5353535353535354,0.5454545454545455,0.5555555555555556,0.5656565656565657,0.5757575757575758,0.5858585858585859,0.595959595959596,0.6060606060606061,0.6161616161616162,0.6262626262626263,0.6363636363636365,0.6464646464646465,0.6565656565656566,0.6666666666666667,0.6767676767676768,0.686868686868687,0.696969696969697,0.7070707070707072,0.7171717171717172,0.7272727272727273,0.7373737373737375,0.7474747474747475,0.7575757575757577,0.7676767676767677,0.7777777777777778,0.787878787878788,0.797979797979798,0.8080808080808082,0.8181818181818182,0.8282828282828284,0.8383838383838385,0.8484848484848485,0.8585858585858587,0.8686868686868687,0.8787878787878789,0.888888888888889,0.8989898989898991,0.9090909090909092,0.9191919191919192,0.9292929292929294,0.9393939393939394,0.9494949494949496,0.9595959595959597,0.9696969696969697,0.9797979797979799,0.98989898989899,1.0],"xaxis":"x","y":[0.0,0.03172793349806765,0.0634239196565645,0.09505604330418266,0.12659245357374926,0.15800139597334992,0.1892512443604102,0.22031053278654067,0.2511479871810792,0.28173255684142967,0.3120334456984871,0.34202014332566877,0.3716624556603275,0.4009305354066137,0.4297949120891717,0.4582265217274104,0.4861967361004687,0.5136773915734064,0.5406408174555976,0.5670598638627707,0.5929079290546405,0.6181589862206052,0.6427876096865394,0.6667690005162916,0.6900790114821119,0.7126941713788629,0.7345917086575333,0.7557495743542583,0.7761464642917569,0.7957618405308321,0.8145759520503357,0.8325698546347714,0.8497254299495144,0.8660254037844387,0.8814533634475821,0.8959937742913359,0.9096319953545183,0.9223542941045815,0.9341478602651067,0.9450008187146685,0.954902241444074,0.963842158559942,0.9718115683235417,0.9788024462147787,0.9848077530122081,0.9898214418809327,0.9938384644612541,0.9968547759519424,0.998867339183008,0.9998741276738751,0.9998741276738751,0.998867339183008,0.9968547759519424,0.9938384644612541,0.9898214418809327,0.984807753012208,0.9788024462147786,0.9718115683235417,0.9638421585599422,0.9549022414440739,0.9450008187146685,0.9341478602651067,0.9223542941045814,0.9096319953545182,0.8959937742913359,0.8814533634475821,0.8660254037844385,0.8497254299495143,0.8325698546347712,0.8145759520503358,0.795761840530832,0.7761464642917568,0.7557495743542583,0.7345917086575331,0.7126941713788627,0.6900790114821119,0.6667690005162917,0.6427876096865395,0.6181589862206051,0.5929079290546404,0.5670598638627704,0.5406408174555974,0.5136773915734063,0.4861967361004687,0.4582265217274105,0.4297949120891714,0.4009305354066136,0.37166245566032713,0.3420201433256685,0.31203344569848696,0.28173255684142967,0.2511479871810793,0.22031053278654036,0.18925124436041008,0.15800139597334945,0.12659245357374938,0.09505604330418288,0.0634239196565644,0.031727933498067656,1.2246467991473532e-16],"yaxis":"y","type":"scatter"},{"legendgroup":"Numerical solution","line":{"color":"blue","dash":"solid"},"marker":{"symbol":"circle"},"mode":"lines","name":"Numerical solution","orientation":"v","showlegend":true,"x":[0.0,0.010101010101010102,0.020202020202020204,0.030303030303030304,0.04040404040404041,0.05050505050505051,0.06060606060606061,0.07070707070707072,0.08080808080808081,0.09090909090909091,0.10101010101010102,0.11111111111111112,0.12121212121212122,0.13131313131313133,0.14141414141414144,0.15151515151515152,0.16161616161616163,0.17171717171717174,0.18181818181818182,0.19191919191919193,0.20202020202020204,0.21212121212121213,0.22222222222222224,0.23232323232323235,0.24242424242424243,0.25252525252525254,0.26262626262626265,0.27272727272727276,0.2828282828282829,0.29292929292929293,0.30303030303030304,0.31313131313131315,0.32323232323232326,0.33333333333333337,0.3434343434343435,0.3535353535353536,0.36363636363636365,0.37373737373737376,0.38383838383838387,0.393939393939394,0.4040404040404041,0.4141414141414142,0.42424242424242425,0.43434343434343436,0.4444444444444445,0.4545454545454546,0.4646464646464647,0.4747474747474748,0.48484848484848486,0.494949494949495,0.5050505050505051,0.5151515151515152,0.5252525252525253,0.5353535353535354,0.5454545454545455,0.5555555555555556,0.5656565656565657,0.5757575757575758,0.5858585858585859,0.595959595959596,0.6060606060606061,0.6161616161616162,0.6262626262626263,0.6363636363636365,0.6464646464646465,0.6565656565656566,0.6666666666666667,0.6767676767676768,0.686868686868687,0.696969696969697,0.7070707070707072,0.7171717171717172,0.7272727272727273,0.7373737373737375,0.7474747474747475,0.7575757575757577,0.7676767676767677,0.7777777777777778,0.787878787878788,0.797979797979798,0.8080808080808082,0.8181818181818182,0.8282828282828284,0.8383838383838385,0.8484848484848485,0.8585858585858587,0.8686868686868687,0.8787878787878789,0.888888888888889,0.8989898989898991,0.9090909090909092,0.9191919191919192,0.9292929292929294,0.9393939393939394,0.9494949494949496,0.9595959595959597,0.9696969696969697,0.9797979797979799,0.98989898989899,1.0],"xaxis":"x","y":[0.0,0.02620845649638624,0.05241691299277244,0.07862536948915867,0.10483382598554487,0.13104228248193112,0.15725073897831726,0.18345919547470355,0.2096676519710897,0.23587610846747598,0.26208456496386223,0.28829302146024843,0.3145014779566346,0.3407099344530209,0.3669183909494071,0.3931268474457933,0.4193353039421795,0.44554376043856575,0.47175221693495195,0.49796067343133815,0.5241691299277245,0.5503775864241106,0.5765860429204969,0.602794499416883,0.6290029559132692,0.6552114124096555,0.6814198689060418,0.7076283254024279,0.7338367818988142,0.7600452383952003,0.7862536948915866,0.8124621513879727,0.838670607884359,0.864879064380745,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807452,0.8648790643807452,0.8648790643807452,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807451,0.8648790643807452,0.8648790643807451,0.8648790643807451,0.8648790643807452,0.8648790643807451,0.8648790643807452,0.8648790643807451,0.8648790643807451,0.8648790643807452,0.8648790643807452,0.8648790643807452,0.8648790643807452,0.8648790643807449,0.8386706078843589,0.8124621513879725,0.7862536948915865,0.7600452383952,0.7338367818988139,0.7076283254024278,0.6814198689060414,0.6552114124096553,0.6290029559132689,0.6027944994168828,0.5765860429204968,0.5503775864241103,0.5241691299277244,0.49796067343133793,0.47175221693495173,0.44554376043856536,0.41933530394217916,0.39312684744579335,0.3669183909494069,0.34070993445302056,0.3145014779566342,0.2882930214602484,0.26208456496386195,0.2358761084674756,0.20966765197108972,0.18345919547470338,0.15725073897831698,0.13104228248193056,0.10483382598554479,0.07862536948915835,0.052416912992772006,0.02620845649638615,0.0],"yaxis":"y","type":"scatter"}],                        {"template":{"data":{"histogram2dcontour":[{"type":"histogram2dcontour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"choropleth":[{"type":"choropleth","colorbar":{"outlinewidth":0,"ticks":""}}],"histogram2d":[{"type":"histogram2d","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmap":[{"type":"heatmap","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmapgl":[{"type":"heatmapgl","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"contourcarpet":[{"type":"contourcarpet","colorbar":{"outlinewidth":0,"ticks":""}}],"contour":[{"type":"contour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"surface":[{"type":"surface","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"mesh3d":[{"type":"mesh3d","colorbar":{"outlinewidth":0,"ticks":""}}],"scatter":[{"fillpattern":{"fillmode":"overlay","size":10,"solidity":0.2},"type":"scatter"}],"parcoords":[{"type":"parcoords","line":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolargl":[{"type":"scatterpolargl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"scattergeo":[{"type":"scattergeo","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolar":[{"type":"scatterpolar","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"scattergl":[{"type":"scattergl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatter3d":[{"type":"scatter3d","line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattermapbox":[{"type":"scattermapbox","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterternary":[{"type":"scatterternary","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattercarpet":[{"type":"scattercarpet","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"type":"carpet"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}],"barpolar":[{"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"pie":[{"automargin":true,"type":"pie"}]},"layout":{"autotypenumbers":"strict","colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"hovermode":"closest","hoverlabel":{"align":"left"},"paper_bgcolor":"white","plot_bgcolor":"#E5ECF6","polar":{"bgcolor":"#E5ECF6","angularaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"radialaxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"ternary":{"bgcolor":"#E5ECF6","aaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"baxis":{"gridcolor":"white","linecolor":"white","ticks":""},"caxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]]},"xaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"yaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"scene":{"xaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"yaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"zaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"geo":{"bgcolor":"white","landcolor":"#E5ECF6","subunitcolor":"white","showland":true,"showlakes":true,"lakecolor":"white"},"title":{"x":0.05},"mapbox":{"style":"light"}}},"xaxis":{"anchor":"y","domain":[0.0,1.0],"title":{"text":"x"}},"yaxis":{"anchor":"x","domain":[0.0,1.0],"title":{"text":"value"}},"legend":{"title":{"text":"variable"},"tracegroupgap":0},"title":{"text":"Manufactured solution of Poisson's equation"},"hovermode":"x","autosize":false,"width":800,"height":500},                        {"responsive": true}                    ).then(function(){

var gd = document.getElementById('9795c319-ae78-4999-8578-d00e1dfa8c61');
var x = new MutationObserver(function (mutations, observer) {{
        var display = window.getComputedStyle(gd).display;
        if (!display || display === 'none') {{
            console.log([gd, 'removed!']);
            Plotly.purge(gd);
            observer.disconnect();
        }}
}});

// Listen for the removal of the full notebook cells
var notebookContainer = gd.closest('#notebook-container');
if (notebookContainer) {{
    x.observe(notebookContainer, {childList: true});
}}

// Listen for the clearing of the current output cell
var outputEl = gd.closest('.output');
if (outputEl) {{
    x.observe(outputEl, {childList: true});
}}

                        })                };                });            </script>        </div>

