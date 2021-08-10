from ctypes import *
from numpy.random import normal
import numpy as np
from dolfin import *
from mesh_generation import sphere_mesh
from utils  import solve_problem
from Problem import Problem 

L = 1
beta = 0.51
kappa =1.0
k = 2
h = 0.01

path_to_c = '../fast_spher_harms.so'
sph = CDLL(path_to_c)
sph.sin_term.restype = c_float
sph.cos_term.restype = c_float
sph.fast_grf.restype = c_float
sph.fast_grf.argtypes = (
    c_int,c_float,c_float,c_float, POINTER(c_float))
def GRF(noterms,rands,x,y,z): #function to evaluate GRF in a point
    return(sph.fast_grf(c_int(L),c_float(x),c_float(y),c_float(z),rands))




def wn(rands,L):
    pf = (c_float*len(rands))(*rands)
    prands = cast(pf,POINTER(c_float))
    class right_hand_side(UserExpression):
        def eval(self,value,x):
            #pf = (c_float*len(rands))(*rands)
            #prands = cast(pf, POINTER(c_float))
            value[0]= GRF(L,prands,x[0],x[1],x[2])
        def value_shape(self):
            return()
    return(right_hand_side(degree=2))


surfaces = [sphere_mesh(h)]
Vs = [FunctionSpace(surfaces[0],FiniteElement("Lagrange",surfaces[0].ufl_cell(),1))]
us = [TrialFunction(Vs[0])]
vs = [TestFunction(Vs[0])]
As = [assemble(inner(grad(us[0]),grad(vs[0]))*dx)]      
Bs = [assemble(us[0]*vs[0]*dx)]

def surf():
    return surfaces[0]
def fenstuff():
    return (Vs[0],us[0],vs[0])
def AB():
    return(As[0],Bs[0])
prob = Problem(surf, fenstuff,wn,beta, k,kappa,AB)
rands = normal(0,1,[L+1,L+1,2]).flatten()
prob.set_rand(rands)



sol = Function(Vs[0])

sol.vector()[:] = solve_problem(prob,5)

file = File('field.pvd')
file << sol
