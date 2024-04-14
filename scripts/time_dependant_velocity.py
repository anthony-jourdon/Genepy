import numpy as np
import bcpy as bp

def main():
  # domain
  O = np.array([0,0],     dtype=np.float64)     # Origin
  L = np.array([600e3,300e3], dtype=np.float64) # Length
  n = np.array([16,16],   dtype=np.int32)       # size
  # cm/a to m/s conversion
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0)
  r_angle = np.deg2rad(0.0)                     # rotation angle
  # velocity phase 1
  u_norm_1  = 1.0 * cma2ms                      # norm of the velocity on boundary
  u_angle_1 = np.deg2rad(90.0)                  # angle of the velocity vector
  u_dir_1   = 0                                 # direction in which velocity varies (0 (x) or 1 (z))
  u_type_1  = "extension"                       # extension or compression (because norm > 0)
  # velocity phase 2
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

  d = bp.Domain(minCoor=O,maxCoor=L,size=n,referential_angle=r_angle)
  bc1 = bp.BoundaryConditions(d,u_norm=u_norm_1,u_angle=u_angle_1,variation_dir=u_dir_1,velocity_type=u_type_1)
  bc2 = bp.BoundaryConditions(d,u_norm=u_norm_2,u_angle=u_angle_2,variation_dir=u_dir_2,velocity_type=u_type_2)
  utils = bp.BCutils(d)

  # step 1: evaluate steady-state velocity during phase 1 & phase 2 (no time involved yet)
  # symbolic evaluation
  uxe,uze = bc1.evaluate_velocity(d.sym_coor[0],d.sym_coor[1])
  uxc,uzc = bc2.evaluate_velocity(d.sym_coor[0],d.sym_coor[1])
  uOe,uLe = bc1.evaluate_boundary_velocity()
  uOc,uLc = bc2.evaluate_boundary_velocity()
  # numerical evaluation
  uxe_n,uze_n = bc1.evaluate_velocity(d.num_coor[0],d.num_coor[1])
  uxc_n,uzc_n = bc2.evaluate_velocity(d.num_coor[0],d.num_coor[1])
  # for 1D plot
  ue = -u_norm_1
  uc = u_norm_2

  # step 2: pass the steady-state velocity to the time evaluator
  # symbolic evaluation
  ux_t,uz_t = utils.symbolic_velocity_inversion((uxe,uze),(uxc,uzc),breakpoints,slopes)
  print('###### uO ######')
  uOx_t,uOz_t = utils.symbolic_velocity_inversion(uOe,uOc,breakpoints,slopes)
  print('###### uL ######')
  uLx_t,uLz_t = utils.symbolic_velocity_inversion(uLe,uLc,breakpoints,slopes)
  # paraview output
  #utils.paraview_velocity_inversion((uxe_n,uze_n),(uxc_n,uzc_n),breakpoints,slopes,time_pv,root,pvd_fname)
  # for 1D plot
  utils.plot_1d_velocity_inversion(ue,uc,breakpoints,slopes,time_plot)

if __name__ == "__main__":
  main()