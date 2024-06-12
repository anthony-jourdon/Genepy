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
  phase_1   = gp.VelocityLinear(Domain,velo1["norm"],velo1["dir"],velo1["type"],velo1["angle"])

  # velocity phase 2 : Compression at 45°
  velo2 = {"norm"  : 3.0 * cma2ms,     # norm of the velocity on boundary
           "angle" : np.deg2rad(45.0), # angle of the velocity vector
           "dir"   : "x",              # direction in which velocity varies 
           "type"  : "compression"     # extension or compression (because norm > 0)
          }
  phase_2   = gp.VelocityLinear(Domain,velo2["norm"],velo2["dir"],velo2["type"],velo2["angle"])

  # time inversion parameters
  Ma2s = (3600.0 * 24.0 * 365.0) * 1e6
  t1 = 5.0 * Ma2s
  t2 = 7.0 * Ma2s
  breakpoints = np.array([ t1, t2 ], dtype=np.float64)  # breakpoints in time (where atan(t-t0)=0 )
  slopes = np.array([ 5e-13, 5e-13 ], dtype=np.float64) # rate of change of the atan() function

  bc_inv = gp.VelocityInversion(Domain,phase_1,phase_2,breakpoints,slopes)
  # space and time dependant velocity function and its gradient
  u,grad_u = bc_inv.evaluate_velocity_and_gradient_symbolic()
  print(bc_inv.report_symbolic_functions(u,grad_u))

  t0 = bc_inv.get_time_zero_velocity(report=True)
  print(f"Time at which velocity is 0: {t0/Ma2s} Myr")

  #t0 = bc_inv.get_time_zero_velocity()
  #print(t0)

  time_1d = np.linspace(0, 20, 201) * Ma2s # time array for plots
  bc_inv.plot_1D_velocity(time_1d)

  """
  time = np.linspace(0, 20, 21) * Ma2s # time array for plots
  root = "./"
  pvd = "timeseries.pvd"
  writer = gp.WriteVTS(Domain)
  bc_inv.paraview_velocity_inversion(writer,time,root,pvd)
  """

if __name__ == "__main__":
  main()