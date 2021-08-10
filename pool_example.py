
from numpy.random import normal
import numpy as np

from dolfin import *
from ctypes import *

from mesh_generation import sphere_mesh
from utils  import solve_problem
from Problem import Problem 

L = 25 # number of KL terms
beta = 0.51 #smoothness parameter 
kappa =1.0 #length scale parameter 
k = 2 #quadrature step size
h = 0.01 #mesh size

path_to_c = './fast_spher_harms.so' #path to right hand side
sph = CDLL(path_to_c)
sph.sin_term.restype = c_float #set types 
sph.cos_term.restype = c_float
sph.fast_grf.restype = c_float
sph.fast_grf.argtypes = ( 
    c_int,c_float,c_float,c_float, POINTER(c_float))
def GRF(noterms,rands,x,y,z): #function to evaluate GRF in a point
    return(sph.fast_grf(c_int(L),c_float(x),c_float(y),c_float(z),rands))




def wn(rands,L): #create right hand side class 
    pf = (c_float*len(rands))(*rands)
    prands = cast(pf,POINTER(c_float)) #set correct type for cfloat. 
    class right_hand_side(UserExpression):
        def eval(self,value,x):
            value[0]= GRF(L,prands,x[0],x[1],x[2])
        def value_shape(self):
            return()
    return(right_hand_side(degree=2))


surfaces = [sphere_mesh(h)] #create mesh, function spaces, trial func, test func, assemble matrices etc. 
Vs = [FunctionSpace(surfaces[0],FiniteElement("Lagrange",surfaces[0].ufl_cell(),1))]
us = [TrialFunction(Vs[0])]
vs = [TestFunction(Vs[0])]
As = [assemble(inner(grad(us[0]),grad(vs[0]))*dx)]      
Bs = [assemble(us[0]*vs[0]*dx)]

def surf(): #put in function so that it can be pickled. 
    return surfaces[0]
def fenstuff():
    return (Vs[0],us[0],vs[0])
def AB():
    return(As[0],Bs[0])
prob = Problem(surf, fenstuff,wn,beta, k,kappa,AB,L=L) #create problem 
rands = normal(0,1,[L+1,L+1,2]).flatten() #set random numbers
prob.set_rand(rands) 



sol = Function(Vs[0]) #placeholder for solution 

sol.vector()[:] = solve_problem(prob,5) #solve!

file = File('field.pvd') #save for plotting. 
file << sol
