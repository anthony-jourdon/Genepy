import numpy as np
import genepy as gp

def main_2d():
  # domain
  O = np.array([0,0],     dtype=np.float64)
  L = np.array([600e3,300e3], dtype=np.float64)
  n = np.array([16,16],   dtype=np.int32)
  # velocity
  r_angle = np.deg2rad(15.0)
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0)
  u_norm  = 1.0 * cma2ms
  u_angle = np.deg2rad(90.0) # u_angle \in [pi/2, pi/2]
  u_dir   = "y"
  u_type  = "compression"

  d  = gp.Domain(2,O,L,n)
  r  = gp.Rotation(2,r_angle)
  bc = gp.Velocity(d,u_norm,u_dir,u_type,u_angle,r)

  u,grad_u = bc.evaluate_velocity_and_gradient_symbolic()

  print("Symbolic velocity function")
  print("ux(x,y) =",u[0,0])
  print("uy(x,y) =",u[0,1])
  print("Derivatives of the velocity function")
  print("dux/dx =",grad_u[0,0])
  print("dux/dy =",grad_u[0,1])
  print("duy/dx =",grad_u[1,0])
  print("duy/dy =",grad_u[1,1])

  u_num = bc.evaluate_velocity_numeric()
  print(u_num[:,2])

  point_data = {"u": u_num}
  w  = gp.WriteVTS(d, vtk_fname="rotated_velocity_2d.vts", point_data=point_data)
  bc.plot_velocity_vts(w)
  #bc.plot_velocity_matplotlib()

def main_3d():
  # domain definition
  O = np.array([0,-250e3,0],    dtype=np.float64) # Origin
  L = np.array([600e3,0,300e3], dtype=np.float64) # Length
  n = np.array([16,16,16],      dtype=np.int32)   # Number of nodes
  # velocity parameters
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0) # cm/a to m/s conversion
  u_norm  = 1.0 * cma2ms                   # horizontal velocity norm
  u_angle = np.deg2rad(90.0)                # velocity angle \in [-pi/2, pi/2]
  u_dir   = "z"                            # direction in which velocity varies
  u_type  = "extension"                    # extension or compression
  # rotation of the referential
  r_angle = np.deg2rad(-15.0)                   # rotation angle
  axis    = np.array([0,1,0], dtype=np.float64) # rotation axis
  # create domain object
  d  = gp.Domain(3,O,L,n)
  # create rotation object
  r  = gp.Rotation(3,r_angle,axis)
  # create boundary conditions object
  bc = gp.Velocity(d,u_norm,u_dir,u_type,u_angle,r)
  
  # evaluate the velocity and its derivatives
  u,grad_u = bc.evaluate_velocity_and_gradient_symbolic()
  uL       = bc.get_velocity_orientation(horizontal=True,normalize=True)
  u_num    = bc.evaluate_velocity_numeric()
  print(bc)
  print(bc.report_symbolic_functions(u,grad_u,uL))
  
  point_data = {"u": u_num}
  w  = gp.WriteVTS(d, vtk_fname="rotated_velocity.vts", point_data=point_data)
  w.write_vts()

  """
  # Sanity check by plugging the numerical coords into the symbolic expression
  u_num = np.zeros(shape=(d.n[0],d.n[1],d.n[2],3), dtype=np.float64)
  u_num[:,:,:,0] = -5.28496533062743e-16*d.num_coor[0] + 1.97237591301416e-15*d.num_coor[2] - 1.37307427033301e-10
  u_num[:,:,:,2] = -1.4161021923681e-16*d.num_coor[0] + 5.28496533062743e-16*d.num_coor[2] - 3.67914141883684e-11
  u_num = np.reshape(u_num,(d.n[0]*d.n[1]*d.n[2],3),order='F')
  """

if __name__ == "__main__":
  #main_2d()
  main_3d()