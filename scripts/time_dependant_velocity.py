import numpy as np
import bcpy as bp

def main():
  # domain
  O = np.array([0,0],     dtype=np.float64)     # Origin
  L = np.array([600e3,300e3], dtype=np.float64) # Length
  n = np.array([16,16],   dtype=np.int32)       # size
  # cm/a to m/s conversion
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0)
  # velocity phase 1
  r_angle_1 = np.deg2rad(0.0)                   # rotation angle
  u_norm_1  = 1.0 * cma2ms                      # norm of the velocity on boundary
  u_angle_1 = np.deg2rad(90.0)                  # angle of the velocity vector
  u_dir_1   = 0                                 # direction in which velocity varies (0 (x) or 1 (z))
  u_type_1  = "extension"                       # extension or compression (because norm > 0)
  # velocity phase 2
  r_angle_2 = np.deg2rad(0.0)
  u_norm_2  = 3.0 * cma2ms
  u_angle_2 = np.deg2rad(45.0)
  u_dir_2   = 0
  u_type_2  = "compression"
  # time inversion parameters
  Ma2s = (3600.0 * 24.0 * 365.0) * 1e6
  t1 = 5.0 * Ma2s
  t2 = 10.0 * Ma2s
  breakpoints = np.array([ t1, t2 ], dtype=np.float64)  # breakpoints in time (where atan(t-t0)=0 )
  slopes = np.array([ 5e-13, 5e-13 ], dtype=np.float64) # rate of change of the atan() function
  time_plot = np.linspace(0, 20, 201) * Ma2s # time array for 1D plot
  time_pv   = np.linspace(0, 20, 21) * Ma2s # time array for paraview output (less step to avoid too much files)
  # paraview output
  root = "./"
  pvd_fname = "timeseries.pvd"

  domain = []
  phase = []
  domain.append(bp.Domain(minCoor=O,maxCoor=L,size=n,referential_angle=r_angle_1))
  domain.append(bp.Domain(minCoor=O,maxCoor=L,size=n,referential_angle=r_angle_2))
  phase.append(bp.Velocity(u_norm_1,u_angle_1,u_dir_1,u_type_1,domain[0]))
  phase.append(bp.Velocity(u_norm_2,u_angle_2,u_dir_2,u_type_2,domain[1]))

  inversion = bp.BCInversion(phase[0],phase[1],breakpoints,slopes)
  writer    = bp.WriteVTS(domain[0])

  inversion.symbolic_velocity_inversion()
  inversion.symbolic_boundary_velocity_inversion()
  
  inversion.paraview_velocity_inversion(writer,time_pv,root,pvd_fname)
  inversion.plot_1d_velocity_inversion(-u_norm_1,u_norm_2,time_plot)


if __name__ == "__main__":
  main()