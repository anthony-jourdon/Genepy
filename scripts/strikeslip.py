import numpy as np
import bcpy as bp

def transpression_7():
  Notche_centre = np.zeros(shape=(2,2), dtype=np.float64)
  # Notch 0 coordinates
  Notche_centre[0,0] = 85900.5  # x
  Notche_centre[0,1] = 123711.5 # z
  # Notch 1 coordinates
  Notche_centre[1,0] = 514099.0 # x
  Notche_centre[1,1] = 176288.0 # z
  return Notche_centre

def transtension(domain,angle,dz):
  domain_centre = 0.5*(domain.O + domain.L)
  wz_centre = np.zeros(shape=(2,2), dtype=np.float64)
  wz_centre[0,1] = domain_centre[1] + dz
  wz_centre[1,1] = domain_centre[1] - dz
  
  wz_centre[0,0] = bp.x_centre_from_angle(wz_centre[0,1],angle,domain_centre)
  wz_centre[1,0] = bp.x_centre_from_angle(wz_centre[1,1],angle,domain_centre)
  return wz_centre

def initial_strain(domain):
  notche_centre = transpression_7()
  #notche_centre = transtension(domain,np.deg2rad(83),25.0e3)
  A = 1.0 #np.random.rand()
  a = 0.5 * 6.0e-5**2
  b = 0.0
  c = 0.5 * 6.0e-5**2
  
  ic = bp.InitialConditions(domain)
  wz_sym,wz_num = ic.evaluate_gaussians(A,a,b,c,notche_centre[:,0],notche_centre[:,1])

  return

def main():
  # domain
  O = np.array([0,0],     dtype=np.float64)
  L = np.array([600e3,300e3], dtype=np.float64)
  n = np.array([64,32],   dtype=np.int32)
  # velocity
  r_angle = np.deg2rad(15.0)
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0)
  u_norm  = 1.0 * cma2ms
  u_angle = np.deg2rad(90.0)
  u_dir   = 1
  u_type  = "extension"

  d = bp.Domain(minCoor=O,maxCoor=L,size=n,referential_angle=r_angle)
  bc0 = bp.BoundaryConditions(d,u_norm=u_norm,u_angle=u_angle,variation_dir=u_dir,velocity_type=u_type)
  bc0.symbolic_derivatives(ccw=True)

  initial_strain(d)

  bc0.plot_velocity(ccw=True)

if __name__ == "__main__":
  main()