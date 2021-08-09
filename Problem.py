from dolfin import *
from numpy import pi,ceil,sin



class Problem:
    def __init__(self,surface,Vuv,rhs,beta,k,kappa,AB,rands=[],L=25):
        self.surface = surface
        self.intbeta = int(beta)
        self.beta = beta-int(beta)
        self.k = k
        self.rhs = rhs
        self.kappa = kappa
        self.rands = rands
        Kmin = int(ceil(pi**2/(4*self.beta*self.k**2)))
        Kplus = int(ceil(pi**2/(4*(1-self.beta)*self.k**2)))
        self.L = L 
        self.Kmin = Kmin
        self.Kplus = Kplus
        self.Vuv = Vuv
        self.AB = AB
    def set_rand(self,rands):
        self.rands = rands
    
    def solve_sequentially(self):
        if(self.rands == []):
            raise(ValueError('Set rands before proceeding'))
        pc = PETScPreconditioner("icc")
        solver = PETScKrylovSolver("cg",pc)
        rhs = self.rhs(self.rands,self.L)
        V,u,v = self.Vuv()
        A, B = self.AB()
        b_temp = assemble(rhs*v*dx)
        w = Function(V)
        if(self.intbeta > 0):
            beta_curr = self.intbeta
            while(beta_curr > 0):
                A_temp = self.kappa**2*B+A
                w = Function(V)
                solver.solve(A_temp,w.vector(),b_temp)
                b_temp = assemble(w*v*dx)
                beta_curr -= 1
            b = b_temp
        else:
            b = b_temp
            
        out = Function(V)
        out_vec = out.vector()[:]
        for l in range(-self.Kmin,self.Kplus+1):
            y = l*self.k
            solver.solve((1+exp(2*y)*self.kappa**2)*B+exp(2*y)*A,out.vector(),b)
            out_vec += exp(2*self.beta*y)*out.vector()[:]
        out_vec *= 2*self.k*sin(pi*self.beta)/pi
        return out_vec 
