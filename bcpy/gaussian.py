import numpy as np
import sympy as sp
from bcpy import utils
from bcpy import domain

import matplotlib.pyplot as plt

class Gaussian(domain.Domain):
  def __init__(self,Domain,ng:int,A,a,b,c,x0,z0) -> None:
    self.ng = ng
    self.A  = np.asarray(A,  dtype=np.float64) # amplitude of the gaussian
    self.a  = np.asarray(a,  dtype=np.float64) # gaussian coefficient
    self.b  = np.asarray(b,  dtype=np.float64) # gaussian coefficient
    self.c  = np.asarray(c,  dtype=np.float64) # gaussian coefficient
    self.x0 = np.asarray(x0, dtype=np.float64) # x coordinate of the gaussian centre
    self.z0 = np.asarray(z0, dtype=np.float64) # z coordinate of the gaussian centre

    # iterate over the attributes and check if they have the same length than ng
    for attribute,value in self.__dict__.items():
      # skip ng attribute
      if attribute == 'ng':
        continue
      # check length
      if value.shape[0] != ng:
        raise RuntimeError(f"Attribute {attribute} must have ng = {ng} elements")

    super(Gaussian, self).__init__(Domain.O,Domain.L,Domain.n,Domain.theta)
    return

  def gaussian_2d(self,A,a,b,c,x,x0,z,z0):
    """
    gaussian_2d(A,a,b,c,x,x0,z,z0)
    Computes a 2D gaussian distribution.

    u = A * exp(-( a*(x-x0)^2 + 2*b*(x-x0)*(z-z0) + c*(z-z0)^2 ) )

    Parameters:
    -----------
    A     : amplitude of the gaussian
    a,b,c : gaussian coefficients
    x0,z0 : center of the gaussian
    x,z   : coordinates

    Returns:
    --------
    u : the field with the 2D gaussian distribution
    """
    if type(x) == sp.core.symbol.Symbol:
      u = A * sp.exp( -( a*(x-x0)*(x-x0) + 2*b*(x-x0)*(z-z0) + c*(z-z0)*(z-z0) ) )
    else:
      u = A * np.exp( -( a*(x-x0)*(x-x0) + 2*b*(x-x0)*(z-z0) + c*(z-z0)*(z-z0) ) )
    return u

  def symbolic_gaussian(self,A,a,b,c,x0,z0):
    g = self.gaussian_2d(A,a,b,c,self.sym_coor[0],x0,self.sym_coor[1],z0)
    return g
  
  def numerical_gaussian(self,A,a,b,c,x0,z0):
    g = self.gaussian_2d(A,a,b,c,self.num_coor[0],x0,self.num_coor[1],z0)
    return g
  
  def evaluate_gaussians(self):
    self.gaussian_sym = np.zeros(shape=(self.ng), dtype=object)
    self.gaussian_num = np.zeros(shape=(self.ng,self.num_coor[0].shape[0],self.num_coor[0].shape[1]), dtype=np.float64)
    for n in range(self.ng):
      print(f'********** Gaussian [{n}] **********')
      self.x0[n],self.z0[n] = self.rotate_referential(self.x0[n],self.z0[n],ccw=True)
      print('Centre:')
      print('[',self.x0[n],self.z0[n],']')
      self.gaussian_sym[n] = self.symbolic_gaussian(self.A[n],self.a[n],self.b[n],self.c[n],self.x0[n],self.z0[n])
      self.gaussian_num[n] = self.numerical_gaussian(self.A[n],self.a[n],self.b[n],self.c[n],self.x0[n],self.z0[n])
      print('Equation:')
      print(self.gaussian_sym[n])
    return 
  
  def plot_gaussians(self):
    _, ax = plt.subplots()
    field = np.zeros(shape=(self.num_coor[0].shape[0],self.num_coor[0].shape[1]), dtype=np.float64)
    for n in range(self.ng):
      field += self.gaussian_num[n,:,:]
    g = ax.contourf(self.num_coor[0],self.num_coor[1],field,100,cmap='magma')
    ax.axis('equal')
    plt.colorbar(g,ax=ax)
    plt.draw()
    return
  
def test():
  # Domain
  O = np.array([ 0.0, 0.0 ],    dtype=np.float64)
  L = np.array([ 600e3, 300e3 ],dtype=np.float64)
  n = np.array([ 64, 32 ],      dtype=np.int32)
  domain_rotation = np.deg2rad(15.0)
  Domain = domain.Domain(O,L,n,domain_rotation)
  # Gaussians
  ng = np.int32(2) # number of gaussians
  A = np.array([1.0, 1.0],dtype=np.float64)
  # shape
  coeff = 0.5 * 6.0e-5**2
  a = np.array([coeff, coeff], dtype=np.float64)
  b = np.array([0.0, 0.0],     dtype=np.float64)
  c = np.array([coeff, coeff], dtype=np.float64)
  # position
  dz    = 25.0e3
  angle = np.deg2rad(83.0)
  domain_centre = 0.5*(Domain.O + Domain.L)
  
  x0 = np.zeros(shape=(ng), dtype=np.float64)
  z0 = np.array([domain_centre[1] - dz, 
                 domain_centre[1] + dz], dtype=np.float64)
  
  x0[0] = utils.x_centre_from_angle(z0[0],angle,domain_centre)
  x0[1] = utils.x_centre_from_angle(z0[1],angle,domain_centre)

  GWZ = Gaussian(Domain,ng,A,a,b,c,x0,z0)
  GWZ.evaluate_gaussians()
  GWZ.plot_gaussians()
  plt.show()
  return

if __name__ == "__main__":
  test()