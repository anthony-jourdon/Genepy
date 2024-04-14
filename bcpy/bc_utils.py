import os
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

class BCutils:
  def __init__(self,Domain) -> None:
    self.domain = Domain
    return
  
  def sum_functions(self,user_func,user_args):
    u = user_func[0](*(user_args[0]))
    nf = len(user_func)
    if nf > 1:
      for n in range(1,nf):
        u += user_func[n](*(user_args[n]))
    return u
  
  def velocity_inversion(self,time,breakpoint,slope,bound):
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
    user_args = [(time,breakpoint[0],slope[0],bound[0]), (time,breakpoint[1],slope[1],bound[1])]
    u_t = self.sum_functions(user_func,user_args)
    return u_t
  
  def symbolic_velocity_inversion(self,u_phase1,u_phase2,breakpoints,slopes):
    time = sp.symbols('t')

    ux_bounds = [ u_phase1[0], u_phase2[0] ]
    uz_bounds = [ u_phase1[1], u_phase2[1] ]
    ux_t = self.velocity_inversion(time,breakpoints,slopes,ux_bounds)
    uz_t = self.velocity_inversion(time,breakpoints,slopes,uz_bounds)
    print('###### velocity function ######')
    print('ux_t =',ux_t)
    print('uz_t =',uz_t)
    duxdx = ux_t.diff(self.domain.sym_coor[0])
    duxdz = ux_t.diff(self.domain.sym_coor[1])
    duzdx = uz_t.diff(self.domain.sym_coor[0])
    duzdz = uz_t.diff(self.domain.sym_coor[1])
    print('###### velocity derivatives ######')
    print('dux/dx =',duxdx)
    print('dux/dz =',duxdz)
    print('duz/dx =',duzdx)
    print('duz/dz =',duzdz)
    return ux_t,uz_t

  def paraview_velocity_inversion(self,u_phase1,u_phase2,breakpoints,slopes,time,root:str,pvd:str):
    u = np.zeros(shape=(self.domain.n[0]*self.domain.n[1],3), dtype=np.float64)

    ux_bounds = [ u_phase1[0], u_phase2[0] ]
    uz_bounds = [ u_phase1[1], u_phase2[1] ]

    pvd_fname = os.path.join(root,pvd)
    step = 0
    self.domain.open_pvd(pvd_fname)
    for t in time:
      output = f"velocity_{step}.vts"

      ux_t = self.velocity_inversion(t,breakpoints,slopes,ux_bounds)
      uz_t = self.velocity_inversion(t,breakpoints,slopes,uz_bounds)
      u[:,0] = ux_t.reshape(self.domain.n[0]*self.domain.n[1],order='F')
      u[:,2] = uz_t.reshape(self.domain.n[0]*self.domain.n[1],order='F')
      point_data = { "u": u }
      self.domain.write_vts(os.path.join(root,output),point_data=point_data)
      self.domain.append_pvd(pvd_fname,t,output)
      step += 1
    self.domain.close_pvd(pvd_fname)
    return
  
  def plot_1d_velocity_inversion(self,u1,u2,breakpoints,slopes,time):
    bounds = [u1, u2]
    u_t = []
    for t in time:
      u_t.append(self.velocity_inversion(t,breakpoints,slopes,bounds))
    u_t = np.asarray(u_t, dtype=np.float32)
    fig,ax = plt.subplots()
    ax.plot(time,u_t)
    plt.show()
    return