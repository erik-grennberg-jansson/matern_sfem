from dolfin import *
from numpy import pi,ceil,sin

#from multiprocessing import Pool
from pathos.multiprocessing import ProcessingPool as Pool
pc = PETScPreconditioner("icc")
solver = PETScKrylovSolver("cg",pc)

               
def solve_subproblem(j,A,B,b,V,k,kappa,beta):
    y = j*k
    w = Function(V)
    solver.solve((1+exp(2*y)*kappa**2)*B+exp(2*y)*A,w.vector(),b) #solve
    return(w.vector()[:]*exp(2*beta*y))


         
def subsolver(args):
    low = int(args[0])
    high = int(args[1])
    prob = args[2]
    surface = prob.surface()
    if(prob.rands == []):
        rhs = prob.rhs()
    else:
        rhs = prob.rhs(prob.rands,prob.L)

    V = prob.Vuv()[0]
    u = prob.Vuv()[1]
    v = prob.Vuv()[2]
    A,B = prob.AB()
    if(prob.intbeta>0):
        a_temp = prob.kappa**2*B+A
        L_temp = assemble(rhs*v*dx)
        w = Function(V)
        solver.solve(a_temp,w,L_temp)
        b = assemble(w*v*dx)
    else:
        b = assemble(rhs*v*dx)
    beta = prob.beta
    kappa = prob.kappa
    k = prob.k
    subsolve = lambda i: solve_subproblem(i,A,B,b,V,k,kappa,beta)
    return(sum(map(subsolve,range(low,high+1))))

        
def solve_problem(problem,procs):
    '''
    if (problem.intbeta>2):
        raise(NotImplementedError('For multipool only beta>2 implemented so far'))
        
    if(problem.rands == []):
        raise(ValueError('Set rands before proceeding'))
    '''
    jobs = []
    sizeSegment = (problem.Kplus+problem.Kmin+1-((problem.Kplus+problem.Kmin+1) % procs))/procs
    for i in range(0, procs):
        jobs.append((i*sizeSegment+1, (i+1)*sizeSegment))
    jobs[-1] = (jobs[-1][0],jobs[-1][1]+((problem.Kplus+problem.Kmin+1)% procs))
    jobs = [(job[0]-(problem.Kmin+1),job[1]-(problem.Kmin+1),problem) for job in jobs ]
    Pool_main = Pool(procs)
    pool = Pool_main.map(subsolver,jobs)
    result = sum(pool)*2*problem.k*sin(pi*problem.beta)/pi
    Pool_main.close()
    return(result)

class normal(UserExpression):
    def eval(self,value,x):
        value[0] = x[0]
        value[1] = x[1]
        value[3] = x[2]
    def value_shape(self):
        return(3,)
