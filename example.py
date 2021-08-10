#necessary imports 
from numpy.random import normal 
import numpy as np
from dolfin import *
from mesh_generation import sphere_mesh
from utils import solve_problem
from field_sfem import problem_const 



L = 100 #number of terms in KL-expansion
beta = 0.51 #smoothness parameter
kappa = 1.0 #length scale parameter
k = 2 #quadrature step size 
h = 0.01 # approx. mesh size 

Problem = problem_const(L,beta,kappa,k,h) #create Problem
 
Problem.set_rand(normal(0,1,[L+1,L+1,2]).flatten()) #set random numbers 

surface = sphere_mesh(h) #create surface 
V= FunctionSpace(surface,FiniteElement("Lagrange",surface.ufl_cell(),1)) #FEniCS function space on surface. 


sol = Function(V) #placeholder for solution. 
print("solving")
sol.vector()[:] = Problem.solve_sequentially() #solve sequentially

print("solved")

print("saving")
file = File('field.pvd') #save for plotting in external software
file << sol 
