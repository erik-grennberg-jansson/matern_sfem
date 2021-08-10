# Matern_sfem 
## Generation of spherical Whittle-Matérn random fields using SFEM 

This package uses surface finite elements combined with a recursion and a Dunford-Taylor quadrature in order to approximate solutions to random partial differential equations of the form 

<a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;(\kappa^2-\Delta_{\mathbb{S}^2})^\beta&space;u&space;=&space;\mathcal{W}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\large&space;(\kappa^2-\Delta_{\mathbb{S}^2})^\beta&space;u&space;=&space;\mathcal{W}" title="\large (\kappa^2-\Delta_{\mathbb{S}^2})^\beta u = \mathcal{W}," /></a>

by means of a surface finite element approximation described in [the paper "Surface finite element approximation of spherical Whittle-Matérn Gaussian random fields"](https://arxiv.org/abs/2102.08822) by Erik Jansson, Mihály Kovács and Annika Lang. 

## Python code 

This Python code is intended to be used with FEniCS 2019.1 and requires an installation of FEniCS 2019.1 to run. Both a naive sequential solver and a solver using the multiprocessing package are included. For information on FEniCS, please see the [homepage of the project](https://fenicsproject.org/). 

Please note that this code is intended to accompany the above paper and as such may not be updated to future releases of FEniCS. 

## Example of usage

For full code, see the file example.py. 

Begin with the necessary imports. Note that while the mesh in this example is created at instantiation, in principle, meshes can also be user-supplied and as such does not need to be created at runtime. 

```python
import numpy as np 
from numpy import normal 
from dolfin import *
from mesh_generation import sphere_mesh 
from utils import solve_problem 
from Problem import Problem 
from field_sfem import problem_const
```
Then decide on parameters. 
```python
L = 25 #number of terms in KL-expansion 
beta = 1.98 #smoothness parameter 
kappa = 1.0 #length scale parameter
k = 2  #quadrature step size.
h = 0.001 #mesh size 
```
If only L is provided, the others will be configured against L as described in the paper. 

Given the parameters, we can create the random field object and set the random numbers. 

```python
Problem = problem_const(L,beta,kappa,k,h)
Problem.set_rand(normal(0,1,[L+1,L+1,2]).flatten())
```
In order to create a FEniCS object which the solution can be saved to, begin by creating a mesh and function space on which a function object is instansiated. 
```python
surface = sphere_mesh(h)
V= FunctionSpace(surface,FiniteElement("Lagrange",surface.ufl_cell(),1))
sol = Function(V)
```
One can then solve and export the solution for plotting in for instance paraview, or further use in a FEniCS-based program. 
```python
sol.vector()[:] = Problem.solve_sequentially()
file = File('field.pvd')
file << sol
```




 
