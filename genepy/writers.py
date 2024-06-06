#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: writers.py
#
#  This file is part of Genepy.
#
#  Genepy is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, either 
#  version 3 of the License, or any later version.
#
#  Genepy is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with Genepy. 
#  If not, see <https://www.gnu.org/licenses/>.
#====================================================================================================


import numpy as np
import pyvista as pvs
from .initial_conditions import domain

class WriteVTS(domain.Domain):
  def __init__(self, Domain, vtk_fname=None, point_data=None, cell_data=None, pvd_fname=None):
    self.vtk_fname  = vtk_fname
    self.point_data = point_data
    self.cell_data  = cell_data
    self.pvd_fname  = pvd_fname
    super(WriteVTS, self).__init__(Domain.dim,Domain.O,Domain.L,Domain.n,coor=Domain.num_coor)

  def write_vts(self):
    coor = list(self.num_coor)
    if self.dim == 2: coor.append(np.zeros(shape=(self.n), dtype=np.float64))
    mesh = pvs.StructuredGrid(*coor)
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
