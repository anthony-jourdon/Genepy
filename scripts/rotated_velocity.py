import numpy as np
import bcpy as bp

def main_2d():
  # domain
  O = np.array([0,0],     dtype=np.float64)
  L = np.array([600e3,300e3], dtype=np.float64)
  n = np.array([16,16],   dtype=np.int32)
  # velocity
  r_angle = np.deg2rad(15.0)
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0)
  u_norm  = 2.0 * cma2ms
  u_angle = np.deg2rad(90.0) # u_angle \in [pi/2, pi/2]
  u_dir   = 1
  u_type  = "compression"

  d  = bp.Domain(2,O,L,n)
  r  = bp.Rotation(2,r_angle)
  bc = bp.BoundaryConditions(u_norm,u_angle,u_dir,u_type,d,r)

  u,grad_u = bc.evaluate_velocity_and_derivatives_symbolic()

  print("Symbolic velocity function")
  print("ux(x,z) =",u[0,0])
  print("uz(x,z) =",u[0,1])
  print("Derivatives of the velocity function")
  print("dux/dx =",grad_u[0,0])
  print("dux/dz =",grad_u[0,1])
  print("duz/dx =",grad_u[1,0])
  print("duz/dz =",grad_u[1,1])

  bc.plot_velocity_matplotlib()

def main_3d():
  # domain
  O = np.array([0,-250e3,0],    dtype=np.float64)
  L = np.array([600e3,0,300e3], dtype=np.float64)
  n = np.array([16,16,16],      dtype=np.int32)
  # velocity
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0)
  u_norm  = 2.0 * cma2ms
  u_angle = np.deg2rad(90.0) # u_angle \in [pi/2, pi/2]
  u_dir   = 2
  u_type  = "extension"
  # rotation
  r_angle = np.deg2rad(-15.0)
  axis    = np.array([0,1,0], dtype=np.float64)

  d  = bp.Domain(3,O,L,n)
  r  = bp.Rotation(3,r_angle,axis)
  bc = bp.BoundaryConditions(u_norm,u_angle,u_dir,u_type,d,r)
  u,grad_u = bc.evaluate_velocity_and_derivatives_symbolic()

  print("Symbolic velocity function")
  print("ux(x,y,z) =",u[0,0])
  print("uy(x,y,z) =",u[0,1])
  print("uz(x,y,z) =",u[0,2])
  print("Derivatives of the velocity function")
  print("dux/dx =",grad_u[0,0])
  print("dux/dy =",grad_u[0,1])
  print("dux/dz =",grad_u[0,2])
  print("duy/dx =",grad_u[1,0])
  print("duy/dy =",grad_u[1,1])
  print("duy/dz =",grad_u[1,2])
  print("duz/dx =",grad_u[2,0])
  print("duz/dy =",grad_u[2,1])
  print("duz/dz =",grad_u[2,2])

  u_num = bc.evaluate_velocity_numeric()  
  point_data = {"u": u_num}
  w  = bp.WriteVTS(d, vtk_fname="rotated_velocity.vts", point_data=point_data)
  bc.plot_velocity_vts(w)

if __name__ == "__main__":
  #main_2d()
  main_3d()