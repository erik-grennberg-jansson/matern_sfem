import meshzoo
import meshio
from dolfin import *
from numpy import exp,ceil,floor

def sphere_mesh(h):
    # https://bitbucket.org/fenics-project/dolfin/issues/845/initialize-mesh-from-vertices
    # Based on code by nschloe/Nico Schlömer
    #approximates a h, NOTE NOT EXACT!
    #TODO: Implement a better solution
    n = int(ceil(h**(-1.0157941567337465)*exp(0.22234672775608288))) 
    
    points,cells = meshzoo.icosa_sphere(n)
    editor = MeshEditor()
    mesh = Mesh()
    editor.open(mesh, "triangle", 2, 3)
    editor.init_vertices(points.shape[0])
    editor.init_cells(cells.shape[0])
    for k, point in enumerate(points):
        editor.add_vertex(k, point)
    for k, cell in enumerate(cells):
        editor.add_cell(k, cell)
    editor.close()        
    return mesh

