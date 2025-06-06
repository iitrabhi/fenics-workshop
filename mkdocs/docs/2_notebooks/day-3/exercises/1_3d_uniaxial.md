The accompanying Jupyter notebook can be obtained here [1_3d_uniaxial](../../../src/day-3/exercises/1_3d_uniaxial.ipynb)



## 3D Uniaxial bar

The objective of this exercise is to implement a 3D finite element analysis (FEA) program using FEniCS. You are required to create a FEniCS code that performs the analysis, solve for the displacements and stresses in the bar, and visualize the results.


Materials:
- Steel material properties:
  - Young's Modulus (E): 200 GPa
  - Poisson's Ratio (ν): 0.3
  - Density (ρ): 7850 kg/m^3

Bar Dimensions:
- Length (Lx): 1 meter
- Ly and Lz: 0.02 meters (20 mm)

Boundary Conditions:
- One end of the bar is fixed (fixed boundary condition).
- The other end is subjected to a uniaxial tensile load:
  - Load (F): 10,000 N

Mesh:
- Use a simple 3D mesh for the bar, consisting of tetrahedral elements.


---
Steps:

1. Import FEniCS and Necessary Libraries
2. Define the Geometry and Mesh
3. Define the Material Properties
4. Define the Boundary Conditions
5. Define the Finite Element Function Space
6. Formulate the Governing Equations.
7.  Solve the System
8.  Calculate Stress
9.  Post-Processing and Visualization:

---


```python
from dolfin import *
```


```python
mesh =  BoxMesh(Point(0,0,0), Point(1,1,1),3,3,3)
```


```python
mesh
```




<!DOCTYPE html>
<html>
  <head>
    <meta http-equiv="content-type" content="text/html;charset=UTF-8" />
    <meta name="generator" content="FEniCS/DOLFIN (http://fenicsproject.org)" />
    <title>FEniCS/DOLFIN X3DOM plot</title>
    <script type="text/javascript" src="https://www.x3dom.org/download/x3dom.js"></script>
    <script type="text/javascript" src="https://code.jquery.com/jquery-3.1.0.js"></script>
    <script type="text/javascript" src="https://rawgit.com/plscott/fenics-x3dom/master/x3dom_support.js"></script>
    <link rel="stylesheet" type="text/css" href="https://www.x3dom.org/download/x3dom.css" />
    <link rel="stylesheet" type="text/css" href="https://rawgit.com/plscott/fenics-x3dom/master/x3dom_support.css" />
  </head>
  <body>
    <x3d showStat="false" xmlns="http://www.web3d.org/specifications/x3d-namespace" width="500.000000px" height="400.000000px">
      <scene>
        <shape>
          <appearance>
            <material diffuseColor="1.000000 1.000000 1.000000" emissiveColor="0.000000 0.000000 0.000000" specularColor="0.000000 0.000000 0.000000" ambientIntensity="0" shininess="0.5" transparency="0"></material>
          </appearance>
          <indexedFaceSet solid="false" colorPerVertex="false" coordIndex="0 5 1 -1 0 17 1 -1 0 5 4 -1 0 20 4 -1 0 17 16 -1 0 20 16 -1 1 6 2 -1 1 18 2 -1 1 6 5 -1 1 18 17 -1 2 7 3 -1 2 19 3 -1 2 7 6 -1 2 19 18 -1 3 21 7 -1 3 21 19 -1 4 9 5 -1 4 9 8 -1 4 22 8 -1 4 22 20 -1 5 10 6 -1 5 10 9 -1 6 11 7 -1 6 11 10 -1 7 23 11 -1 7 23 21 -1 8 13 9 -1 8 13 12 -1 8 24 12 -1 8 24 22 -1 9 14 10 -1 9 14 13 -1 10 15 11 -1 10 15 14 -1 11 27 15 -1 11 27 23 -1 12 25 13 -1 12 25 24 -1 13 26 14 -1 13 26 25 -1 14 27 15 -1 14 27 26 -1 16 29 17 -1 16 32 20 -1 16 29 28 -1 16 32 28 -1 17 30 18 -1 17 30 29 -1 18 31 19 -1 18 31 30 -1 19 33 21 -1 19 33 31 -1 20 34 22 -1 20 34 32 -1 21 35 23 -1 21 35 33 -1 22 36 24 -1 22 36 34 -1 23 39 27 -1 23 39 35 -1 24 37 25 -1 24 37 36 -1 25 38 26 -1 25 38 37 -1 26 39 27 -1 26 39 38 -1 28 41 29 -1 28 44 32 -1 28 41 40 -1 28 44 40 -1 29 42 30 -1 29 42 41 -1 30 43 31 -1 30 43 42 -1 31 47 33 -1 31 47 43 -1 32 48 34 -1 32 48 44 -1 33 51 35 -1 33 51 47 -1 34 52 36 -1 34 52 48 -1 35 55 39 -1 35 55 51 -1 36 53 37 -1 36 53 52 -1 37 54 38 -1 37 54 53 -1 38 55 39 -1 38 55 54 -1 40 45 41 -1 40 45 44 -1 41 46 42 -1 41 46 45 -1 42 47 43 -1 42 47 46 -1 44 49 45 -1 44 49 48 -1 45 50 46 -1 45 50 49 -1 46 51 47 -1 46 51 50 -1 48 53 49 -1 48 53 52 -1 49 54 50 -1 49 54 53 -1 50 55 51 -1 50 55 54 -1 ">
            <coordinate point="0 0 0 0.333333 0 0 0.666667 0 0 1 0 0 0 0.333333 0 0.333333 0.333333 0 0.666667 0.333333 0 1 0.333333 0 0 0.666667 0 0.333333 0.666667 0 0.666667 0.666667 0 1 0.666667 0 0 1 0 0.333333 1 0 0.666667 1 0 1 1 0 0 0 0.333333 0.333333 0 0.333333 0.666667 0 0.333333 1 0 0.333333 0 0.333333 0.333333 1 0.333333 0.333333 0 0.666667 0.333333 1 0.666667 0.333333 0 1 0.333333 0.333333 1 0.333333 0.666667 1 0.333333 1 1 0.333333 0 0 0.666667 0.333333 0 0.666667 0.666667 0 0.666667 1 0 0.666667 0 0.333333 0.666667 1 0.333333 0.666667 0 0.666667 0.666667 1 0.666667 0.666667 0 1 0.666667 0.333333 1 0.666667 0.666667 1 0.666667 1 1 0.666667 0 0 1 0.333333 0 1 0.666667 0 1 1 0 1 0 0.333333 1 0.333333 0.333333 1 0.666667 0.333333 1 1 0.333333 1 0 0.666667 1 0.333333 0.666667 1 0.666667 0.666667 1 1 0.666667 1 0 1 1 0.333333 1 1 0.666667 1 1 1 1 1 "></coordinate>
          </indexedFaceSet>
        </shape>
        <shape>
          <appearance>
            <material emissiveColor="0.000000 0.000000 0.000000" specularColor="0.000000 0.000000 0.000000" ambientIntensity="0" shininess="0.5" transparency="0"></material>
          </appearance>
          <indexedLineSet solid="false" colorPerVertex="false" coordIndex="0 1 -1 0 4 -1 0 5 -1 0 16 -1 0 17 -1 0 20 -1 1 2 -1 1 5 -1 1 6 -1 1 17 -1 1 18 -1 2 3 -1 2 6 -1 2 7 -1 2 18 -1 2 19 -1 3 7 -1 3 19 -1 3 21 -1 4 5 -1 4 8 -1 4 9 -1 4 20 -1 4 22 -1 5 6 -1 5 9 -1 5 10 -1 6 7 -1 6 10 -1 6 11 -1 7 11 -1 7 21 -1 7 23 -1 8 9 -1 8 12 -1 8 13 -1 8 22 -1 8 24 -1 9 10 -1 9 13 -1 9 14 -1 10 11 -1 10 14 -1 10 15 -1 11 15 -1 11 23 -1 11 27 -1 12 13 -1 12 24 -1 12 25 -1 13 14 -1 13 25 -1 13 26 -1 14 15 -1 14 26 -1 14 27 -1 15 27 -1 16 17 -1 16 20 -1 16 28 -1 16 29 -1 16 32 -1 17 18 -1 17 29 -1 17 30 -1 18 19 -1 18 30 -1 18 31 -1 19 21 -1 19 31 -1 19 33 -1 20 22 -1 20 32 -1 20 34 -1 21 23 -1 21 33 -1 21 35 -1 22 24 -1 22 34 -1 22 36 -1 23 27 -1 23 35 -1 23 39 -1 24 25 -1 24 36 -1 24 37 -1 25 26 -1 25 37 -1 25 38 -1 26 27 -1 26 38 -1 26 39 -1 27 39 -1 28 29 -1 28 32 -1 28 40 -1 28 41 -1 28 44 -1 29 30 -1 29 41 -1 29 42 -1 30 31 -1 30 42 -1 30 43 -1 31 33 -1 31 43 -1 31 47 -1 32 34 -1 32 44 -1 32 48 -1 33 35 -1 33 47 -1 33 51 -1 34 36 -1 34 48 -1 34 52 -1 35 39 -1 35 51 -1 35 55 -1 36 37 -1 36 52 -1 36 53 -1 37 38 -1 37 53 -1 37 54 -1 38 39 -1 38 54 -1 38 55 -1 39 55 -1 40 41 -1 40 44 -1 40 45 -1 41 42 -1 41 45 -1 41 46 -1 42 43 -1 42 46 -1 42 47 -1 43 47 -1 44 45 -1 44 48 -1 44 49 -1 45 46 -1 45 49 -1 45 50 -1 46 47 -1 46 50 -1 46 51 -1 47 51 -1 48 49 -1 48 52 -1 48 53 -1 49 50 -1 49 53 -1 49 54 -1 50 51 -1 50 54 -1 50 55 -1 51 55 -1 52 53 -1 53 54 -1 54 55 -1 ">
            <coordinate point="0 0 0 0.333333 0 0 0.666667 0 0 1 0 0 0 0.333333 0 0.333333 0.333333 0 0.666667 0.333333 0 1 0.333333 0 0 0.666667 0 0.333333 0.666667 0 0.666667 0.666667 0 1 0.666667 0 0 1 0 0.333333 1 0 0.666667 1 0 1 1 0 0 0 0.333333 0.333333 0 0.333333 0.666667 0 0.333333 1 0 0.333333 0 0.333333 0.333333 1 0.333333 0.333333 0 0.666667 0.333333 1 0.666667 0.333333 0 1 0.333333 0.333333 1 0.333333 0.666667 1 0.333333 1 1 0.333333 0 0 0.666667 0.333333 0 0.666667 0.666667 0 0.666667 1 0 0.666667 0 0.333333 0.666667 1 0.333333 0.666667 0 0.666667 0.666667 1 0.666667 0.666667 0 1 0.666667 0.333333 1 0.666667 0.666667 1 0.666667 1 1 0.666667 0 0 1 0.333333 0 1 0.666667 0 1 1 0 1 0 0.333333 1 0.333333 0.333333 1 0.666667 0.333333 1 1 0.333333 1 0 0.666667 1 0.333333 0.666667 1 0.666667 0.666667 1 1 0.666667 1 0 1 1 0.333333 1 1 0.666667 1 1 1 1 1 "></coordinate>
          </indexedLineSet>
        </shape>
        <viewpoint id="default" position="2.267767 2.267767 2.267767" orientation="-0.7071067812 0.7071067812 0 1" fieldOfView="0.785398" centerOfRotation="0.500000 0.500000 0.500000" zNear="-1" zFar="-1"></viewpoint>
        <viewpoint id="top" position="0.500000 3.000000 0.500000" orientation="-1 0 0 1.5707963267948" fieldOfView="0.785398" centerOfRotation="0.500000 0.500000 0.500000" zNear="-1" zFar="-1"></viewpoint>
        <viewpoint id="bottom" position="0.500000 -2.000000 0.500000" orientation="1 0 0 1.5707963267948" fieldOfView="0.785398" centerOfRotation="0.500000 0.500000 0.500000" zNear="-1" zFar="-1"></viewpoint>
        <viewpoint id="left" position="3.000000 0.500000 0.500000" orientation="0 1 0 1.5707963267948" fieldOfView="0.785398" centerOfRotation="0.500000 0.500000 0.500000" zNear="-1" zFar="-1"></viewpoint>
        <viewpoint id="right" position="-2.000000 0.500000 0.500000" orientation="0 -1 0 1.5707963267948" fieldOfView="0.785398" centerOfRotation="0.500000 0.500000 0.500000" zNear="-1" zFar="-1"></viewpoint>
        <viewpoint id="back" position="0.500000 0.500000 -2.500000" orientation="0 1 0 3.1415926535898" fieldOfView="0.785398" centerOfRotation="0.500000 0.500000 0.500000" zNear="-1" zFar="-1"></viewpoint>
        <viewpoint id="front" position="0.500000 0.500000 3.000000" orientation="0 0 0 1" fieldOfView="0.785398" centerOfRotation="0.500000 0.500000 0.500000" zNear="-1" zFar="-1"></viewpoint>
        <background skyColor="0.950000 0.950000 0.950000"></background>
        <directionalLight ambientIntensity="0" intensity="1"></directionalLight>
      </scene>
    </x3d>
    <div id="menu"><form id="menu-items"><input type="radio" id="button-options" name="menu" checked="" />
        <label for="button-options">Options</label>
        <input type="radio" id="button-summary" name="menu" />
        <label for="button-summary" style="display: none;">Summary</label>
        <input type="radio" id="button-color" name="menu" />
        <label for="button-color" style="display: none;">Color</label>
        <input type="radio" id="button-warp" name="menu" />
        <label for="button-warp" style="display: none;">Warp</label>
        <input type="radio" id="button-viewpoints" name="menu" />
        <label for="button-viewpoints" style="display: none;">Viewpoints</label>
      </form>
      <div id="menu-content"><div id="content-options" for="button-options"><span>Menu Options</span>
          <br />
          <form class="options"><input type="checkbox" id="select-summary" />
            <label for="select-summary">Summary</label>
          </form>
          <br />
          <form class="options"><input type="checkbox" id="select-color" />
            <label for="select-color">Color</label>
          </form>
          <br />
          <form class="options"><input type="checkbox" id="select-warp" />
            <label for="select-warp">Warp</label>
          </form>
          <br />
          <form class="options"><input type="checkbox" id="select-viewpoints" />
            <label for="select-viewpoints">Viewpoints</label>
          </form>
        </div>
        <div id="content-summary" for="button-summary" hidden="">Number of vertices: 64<br />Number of cells: 162</div>
        <div id="content-color" for="button-color" hidden=""><form>
            <input id="color-checkbox" type="checkbox" checked="" />
            <label for="color-checkbox">Show Color</label>
          </form>Current Color Map:<br />
          <span id="min-color-value"></span>
          <span id="color-map"></span>
          <span id="max-color-value"></span>
        </div>
        <div id="content-warp" for="button-warp" hidden=""><form>
            <input id="warp-checkbox" type="checkbox" />
            <label for="warp-checkbox">Warp By Scalar</label>
            <br />
            <input id="warp-slider" type="range" min="0" max="5" step="0.01" value="1" disabled="" />
            <br />
            <label id="warp-slider-val" for="warp-slider">1</label>
          </form>
        </div>
        <div id="content-viewpoints" for="button-viewpoints" hidden=""><span>Viewpoint Options</span>
          <br />
          <button class="viewpoint">front</button>
          <button class="viewpoint">back</button>
          <button class="viewpoint">left</button>
          <br />
          <button class="viewpoint">right</button>
          <button class="viewpoint">top</button>
          <button class="viewpoint">bottom</button>
        </div>
      </div>
    </div>
  </body>
</html>





```python

```
