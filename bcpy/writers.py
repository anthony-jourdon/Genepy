
import numpy as np
import pyvista as pvs
from . import domain

class WriteVTS(domain.Domain):
  def __init__(self, Domain, vtk_fname=None, point_data=None, cell_data=None, pvd_fname=None):
    self.vtk_fname  = vtk_fname
    self.point_data = point_data
    self.cell_data  = cell_data
    self.pvd_fname  = pvd_fname
    super(WriteVTS, self).__init__(Domain.dim,Domain.O,Domain.L,Domain.n)

  def write_vts(self):
    if self.dim == 2:
      y = np.zeros(shape=(self.num_coor[0].shape), dtype=np.float64)
      mesh = pvs.StructuredGrid(self.num_coor[0],y,self.num_coor[1])
    elif self.dim == 3:
      mesh = pvs.StructuredGrid(self.num_coor[0],self.num_coor[1],self.num_coor[2])

    if self.point_data:
      for data in self.point_data:
        mesh.point_data[data] = self.point_data[data]
    if self.cell_data:
      for data in self.cell_data:
        mesh.cell_data[data] = self.cell_data[data]
    #print(mesh)
    mesh.save(self.vtk_fname)
    return
  
  def open_pvd(self):
    with open(self.pvd_fname,"w") as f:
      f.write('<?xml version="1.0"?>\n')
      f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
      f.write('<Collection>\n')
    return 
  
  def close_pvd(self):
    with open(self.pvd_fname,"a") as f:
      f.write('</Collection>\n')
      f.write('</VTKFile>')
    return

  def append_pvd(self,time,vtk_file):
    with open(self.pvd_fname,"a") as f:
      f.write(f'  <DataSet timestep="{time}" file="{vtk_file}"/>\n')
    return