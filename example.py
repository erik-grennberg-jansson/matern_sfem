from numpy.random import normal
import numpy as np
from dolfin import *
from mesh_generation import sphere_mesh
from utils import solve_problem
#from Problem import Problem
from field_sfem import problem_const



L = 100
beta = 0.51
kappa = 1.0
k = 2
h = 0.01

Problem = problem_const(L,beta,kappa,k,h)

Problem.set_rand(normal(0,1,[L+1,L+1,2]).flatten())
print(Problem)
surface = sphere_mesh(h)
V= FunctionSpace(surface,FiniteElement("Lagrange",surface.ufl_cell(),1))


sol = Function(V)
print("solving")
#sol.vector()[:] = Problem.solve_sequentially() #SEQ SOLVE IS SLOW! 
sol.vector()[:] = solve_problem(Problem,5)
print("solved")

print("saving")
file = File('field.pvd')
file << sol
