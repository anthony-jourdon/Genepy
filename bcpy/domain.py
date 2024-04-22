
import numpy as np
import sympy as sp

class Domain:
  def __init__(self,minCoor=np.zeros(shape=(2), dtype=np.float64),
                    maxCoor=np.ones(shape=(2), dtype=np.float64),
                    size=np.array([2, 2], dtype=np.int32),
                    referential_angle=0.0) -> None:
    self.O     = minCoor
    self.L     = maxCoor
    self.n     = size
    self.theta = referential_angle

    self.num_coor = None
    self.sym_coor = None

    self.numerical_coordinates()
    self.symbolic_coordinates()

  def numerical_coordinates(self):
    coor = list()
    for i in range(self.n.shape[0]):
      coor.append( np.linspace(self.O[i],self.L[i],self.n[i], dtype=np.float64) )
    X,Z = np.meshgrid(coor[0],coor[1])
    self.num_coor = (X,Z)
    return

  def symbolic_coordinates(self):
    self.sym_coor = sp.symbols('x z')
    return
  
  def rotate_referential(self,x,z,ccw=False):
    # translate
    xT = x - 0.5*(self.L[0] + self.O[0])
    zT = z - 0.5*(self.L[1] + self.O[1])
    # rotate
    if ccw:
      xTR = xT * np.cos(self.theta) - zT * np.sin(self.theta)
      zTR = zT * np.cos(self.theta) + xT * np.sin(self.theta)
    else:
      xTR = xT * np.cos(self.theta) + zT * np.sin(self.theta)
      zTR = zT * np.cos(self.theta) - xT * np.sin(self.theta)
    # translate back
    xR = xTR + 0.5*(self.L[0] + self.O[0])
    zR = zTR + 0.5*(self.L[1] + self.O[1])
    return xR,zR
  
  def write_vts(self,fname,point_data=None,cell_data=None):
    import pyvista as pvs

    y = np.zeros(shape=(self.num_coor[0].shape), dtype=np.float64)
    mesh = pvs.StructuredGrid(self.num_coor[0],y,self.num_coor[1])

    if point_data:
      for data in point_data:
        mesh.point_data[data] = point_data[data]
    if cell_data:
      for data in cell_data:
        mesh.cell_data[data] = cell_data[data]
    #print(mesh)
    mesh.save(fname)
    return
  
  def open_pvd(self,fname):
    with open(fname,"w") as f:
      f.write('<?xml version="1.0"?>\n')
      f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
      f.write('<Collection>\n')
    return 
  
  def close_pvd(self,fname):
    with open(fname,"a") as f:
      f.write('</Collection>\n')
      f.write('</VTKFile>')
    return

  def append_pvd(self,fname,time,vts_file):
    with open(fname,"a") as f:
      f.write(f'  <DataSet timestep="{time}" file="{vts_file}"/>\n')
    return