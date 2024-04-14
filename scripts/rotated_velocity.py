import numpy as np
import bcpy as bp

def main():
  # domain
  O = np.array([0,0],     dtype=np.float64)
  L = np.array([600e3,300e3], dtype=np.float64)
  n = np.array([16,16],   dtype=np.int32)
  # velocity
  r_angle = np.deg2rad(-15.0)
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0)
  u_norm  = 2.0 * cma2ms
  u_angle = np.deg2rad(90.0)
  u_dir   = 1
  u_type  = "compression"

  d = bp.Domain(minCoor=O,maxCoor=L,size=n,referential_angle=r_angle)
  bc0 = bp.BoundaryConditions(d,u_norm=u_norm,u_angle=u_angle,variation_dir=u_dir,velocity_type=u_type)
  bc0.symbolic_derivatives()
  bc0.plot_velocity()

if __name__ == "__main__":
  main()