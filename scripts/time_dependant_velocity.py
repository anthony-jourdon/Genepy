import numpy as np
import genepy as gp

def main():
  # 3D domain
  dimensions = 3
  O = np.array([0,-100e3,0],    dtype=np.float64) # Origin
  L = np.array([600e3,0,300e3], dtype=np.float64) # Length
  n = np.array([9,9,9],         dtype=np.int32)   # Number of Q1 nodes i.e. elements + 1
  # Create Domain class instance
  Domain = gp.Domain(dimensions,O,L,n)

  # cm/a to m/s conversion
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0)

  # velocity phase 1 : Orthogonal extension
  velo1 = {"norm"  : 1.0 * cma2ms,     # norm of the velocity on boundary
           "angle" : np.deg2rad(90.0), # angle of the velocity vector
           "dir"   : "x",              # direction in which velocity varies 
           "type"  : "extension"       # extension or compression (because norm > 0)
          }
  phase_1   = gp.Velocity(Domain,velo1["norm"],velo1["dir"],velo1["type"],velo1["angle"])

  # velocity phase 2 : Compression at 45°
  velo2 = {"norm"  : 3.0 * cma2ms,     # norm of the velocity on boundary
           "angle" : np.deg2rad(45.0), # angle of the velocity vector
           "dir"   : "x",              # direction in which velocity varies 
           "type"  : "compression"     # extension or compression (because norm > 0)
          }
  phase_2   = gp.Velocity(Domain,velo2["norm"],velo2["dir"],velo2["type"],velo2["angle"])

  # time inversion parameters
  Ma2s = (3600.0 * 24.0 * 365.0) * 1e6
  t1 = 5.0 * Ma2s
  t2 = 10.0 * Ma2s
  breakpoints = np.array([ t1, t2 ], dtype=np.float64)  # breakpoints in time (where atan(t-t0)=0 )
  slopes = np.array([ 5e-13, 5e-13 ], dtype=np.float64) # rate of change of the atan() function

  bc_inv = gp.VelocityInversion(phase_1,phase_2,breakpoints,slopes)
  # space and time dependant velocity function
  u_t = bc_inv.symbolic_velocity_inversion()
  # gradient
  


  time = np.linspace(0, 20, 21) * Ma2s # time array for plots
  root = "./"
  pvd = "timeseries.pvd"
  writer = gp.WriteVTS(Domain)
  bc_inv.paraview_velocity_inversion(writer,time,root,pvd)

  """
  time_plot = np.linspace(0, 20, 201) * Ma2s # time array for 1D plot
  time_pv   = np.linspace(0, 20, 21) * Ma2s # time array for paraview output (less step to avoid too much files)
  # paraview output
  root = "./"
  pvd_fname = "timeseries.pvd"

  domain = []
  phase = []
  domain.append(gp.Domain(minCoor=O,maxCoor=L,size=n,referential_angle=r_angle_1))
  domain.append(gp.Domain(minCoor=O,maxCoor=L,size=n,referential_angle=r_angle_2))
  phase.append(gp.Velocity(u_norm_1,u_angle_1,u_dir_1,u_type_1,domain[0]))
  phase.append(gp.Velocity(u_norm_2,u_angle_2,u_dir_2,u_type_2,domain[1]))

  inversion = gp.BCInversion(phase[0],phase[1],breakpoints,slopes)
  writer    = gp.WriteVTS(domain[0])

  inversion.symbolic_velocity_inversion()
  inversion.symbolic_boundary_velocity_inversion()
  
  inversion.paraview_velocity_inversion(writer,time_pv,root,pvd_fname)
  inversion.plot_1d_velocity_inversion(-u_norm_1,u_norm_2,time_plot)
  """


if __name__ == "__main__":
  main()