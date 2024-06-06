#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  bc-pre-processing
#  filename: bc_inversion.py
#
#  This file is part of bc-pre-processing.
#
#  bc-pre-processing is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, either 
#  version 3 of the License, or any later version.
#
#  bc-pre-processing is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with bc-pre-processing. 
#  If not, see <https://www.gnu.org/licenses/>.
#====================================================================================================

import os
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

class BCInversion:
  def __init__(self,phase1,phase2,breakpoints,slopes) -> None:
    self.phase1      = phase1
    self.phase2      = phase2
    self.breakpoints = breakpoints
    self.slopes      = slopes
    self.time_sym    = sp.symbols('t')
  
  def sum_functions(self,user_func,user_args):
    u = user_func[0](*(user_args[0]))
    nf = len(user_func)
    if nf > 1:
      for n in range(1,nf):
        u += user_func[n](*(user_args[n]))
    return u
  
  def velocity_inversion(self,time,bound):
    """
    velocity_inversion()
    Provides a velocity inversion function by performing the addition of 
    2 arctangent functions.
    
    The solution provided in this function is manufactured but it can be used as an exemple
    to build more complicated cases with more functions
    """
    def phase_1(t,t0,s,b):
      fac = 2.0/np.pi
      func = 0.5*( -fac*b*sp.atan(s*(t-t0)) + b )
      return func
    
    def phase_2(t,t0,s,b):
      fac = 2.0/np.pi
      func = 0.5*( fac*b*sp.atan(s*(t-t0)) + b )
      return func
    
    user_func = [phase_1,phase_2]
    user_args = [(time,self.breakpoints[0],self.slopes[0],bound[0]), (time,self.breakpoints[1],self.slopes[1],bound[1])]
    u_t = self.sum_functions(user_func,user_args)
    return u_t
  
  def symbolic_velocity_inversion(self):
    ux1,uz1 = self.phase1.evaluate_velocity(self.phase1.sym_coor[0],self.phase1.sym_coor[1])
    ux2,uz2 = self.phase2.evaluate_velocity(self.phase2.sym_coor[0],self.phase2.sym_coor[1])

    ux_bounds = [ ux1, ux2 ]
    uz_bounds = [ uz1, uz2 ]

    ux_t = self.velocity_inversion(self.time_sym,ux_bounds)
    uz_t = self.velocity_inversion(self.time_sym,uz_bounds)

    print('###### velocity function ######')
    print('ux_t =',ux_t)
    print('uz_t =',uz_t)
    self.phase1.symbolic_derivatives(ux_t,uz_t)
    return

  def symbolic_boundary_velocity_inversion(self):
    uO1,uL1 = self.phase1.evaluate_boundary_velocity()
    uO2,uL2 = self.phase2.evaluate_boundary_velocity()

    print('###### uO ######')
    uOx_bounds = [ uO1[0], uO2[0] ]
    uOz_bounds = [ uO1[1], uO2[1] ]
    uOx_t = self.velocity_inversion(self.time_sym,uOx_bounds)
    uOz_t = self.velocity_inversion(self.time_sym,uOz_bounds)
    print('x component =', uOx_t)
    print('z component =', uOz_t)
    
    print('###### uL ######')
    uLx_bounds = [ uL1[0], uL2[0] ]
    uLz_bounds = [ uL1[1], uL2[1] ]
    uLx_t = self.velocity_inversion(self.time_sym,uLx_bounds)
    uLz_t = self.velocity_inversion(self.time_sym,uLz_bounds)
    print('x component =', uLx_t)
    print('z component =', uLz_t)
    return
  
  def paraview_velocity_inversion(self,writer,time,root:str,pvd:str):

    u = np.zeros(shape=(writer.n[0]*writer.n[1],3), dtype=np.float64)

    ux1,uz1 = self.phase1.evaluate_velocity(self.phase1.num_coor[0],self.phase1.num_coor[1])
    ux2,uz2 = self.phase2.evaluate_velocity(self.phase2.num_coor[0],self.phase2.num_coor[1])

    ux_bounds = [ ux1, ux2 ]
    uz_bounds = [ uz1, uz2 ]

    writer.pvd_fname = os.path.join(root,pvd)

    step = 0
    writer.open_pvd()
    for t in time:
      output = f"velocity_{step}.vts"
      writer.vtk_fname = os.path.join(root,output)

      ux_t = self.velocity_inversion(t,ux_bounds)
      uz_t = self.velocity_inversion(t,uz_bounds)
      u[:,0] = ux_t.reshape(writer.n[0]*writer.n[1],order='F')
      u[:,2] = uz_t.reshape(writer.n[0]*writer.n[1],order='F')
      writer.point_data = { "u": u }
      writer.write_vts()
      writer.append_pvd(t,output)
      step += 1
    writer.close_pvd()                                                            
    return
  
  def plot_1d_velocity_inversion(self,u1,u2,time):
    bounds = [u1, u2]
    u_t = []
    for t in time:
      u_t.append(self.velocity_inversion(t,bounds))
    u_t = np.asarray(u_t, dtype=np.float32)
    fig,ax = plt.subplots()
    ax.plot(time,u_t)
    plt.show()
    return
