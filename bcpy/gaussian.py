import numpy as np
import sympy as sp
from bcpy import domain
from bcpy import rotation

import matplotlib.pyplot as plt

class Gaussian(domain.Domain,rotation.Rotation):
  def __init__(self,Domain,Rotation,ng:int,A,a,b,c,x0,z0) -> None:
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

    domain.Domain.__init__(self,Domain.dim,Domain.O,Domain.L,Domain.n)
    rotation.Rotation.__init__(self,Domain.dim,Rotation.theta,Rotation.axis)
    return
  
  def __str__(self) -> str:
    s = f'Gaussian distribution for {self.ng} gaussians\n'
    for n in range(self.ng):
      s += f'Gaussian [{n}]\n'
      s += f'\tAmplitude: {self.A[n]}\n'
      s += f'\ta: {self.a[n]}\n'
      s += f'\tb: {self.b[n]}\n'
      s += f'\tc: {self.c[n]}\n'
      s += f'\tCentre: [ {self.x0[n]},{self.z0[n]} ]\n'
    return s

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
    if self.dim == 2:
      x = self.sym_coor[0]
      z = self.sym_coor[1]
    elif self.dim == 3:
      x = self.sym_coor[0]
      z = self.sym_coor[2]
    g = self.gaussian_2d(A,a,b,c,x,x0,z,z0)
    return g
  
  def numerical_gaussian(self,A,a,b,c,x0,z0):
    if self.dim == 2:
      x = self.num_coor[0]
      z = self.num_coor[1]
    elif self.dim == 3:
      x = self.num_coor[0]
      z = self.num_coor[2]
    g = self.gaussian_2d(A,a,b,c,x,x0,z,z0)
    return g
  
  def evaluate_gaussians(self):
    self.gaussian_sym = np.zeros(shape=(self.ng), dtype=object)
    if self.dim == 2:
      self.gaussian_num = np.zeros(shape=(self.ng,self.num_coor[0].shape[0],self.num_coor[0].shape[1]), dtype=np.float64)
    elif self.dim == 3:
      self.gaussian_num = np.zeros(shape=(self.ng,self.num_coor[0].shape[0],self.num_coor[0].shape[2]), dtype=np.float64)
    
    for n in range(self.ng):
      print(f'********** Gaussian [{n}] **********')
      g_centre = np.array([self.x0[n],self.z0[n]],dtype=np.float64).T
      g_centre = self.rotate_referential(g_centre,self.O,self.L)
      self.x0[n] = g_centre[0]
      self.z0[n] = g_centre[1]
      print('Centre:')
      print(f'[ {self.x0[n]},{self.z0[n]} ]')
      self.gaussian_sym[n] = self.symbolic_gaussian(self.A[n],self.a[n],self.b[n],self.c[n],self.x0[n],self.z0[n])
      self.gaussian_num[n] = self.numerical_gaussian(self.A[n],self.a[n],self.b[n],self.c[n],self.x0[n],self.z0[n])
      print('Equation:')
      print(self.gaussian_sym[n])
    return 
  
  def plot_gaussians(self):
    _, ax = plt.subplots()
    if self.dim == 2:
      field = np.zeros(shape=(self.num_coor[0].shape[0],self.num_coor[0].shape[1]), dtype=np.float64)
    elif self.dim == 3:
      field = np.zeros(shape=(self.num_coor[0].shape[0],self.num_coor[0].shape[2]), dtype=np.float64)

    for n in range(self.ng):
      field += self.gaussian_num[n,:,:]
    g = ax.contourf(self.num_coor[0],self.num_coor[1],field,100,cmap='magma')
    ax.axis('equal')
    plt.colorbar(g,ax=ax)
    plt.draw()
    return
  
def test():
  from bcpy import utils
  from bcpy import rotation
  # Domain
  O = np.array([ 0.0, 0.0 ],    dtype=np.float64)
  L = np.array([ 600e3, 300e3 ],dtype=np.float64)
  n = np.array([ 64, 32 ],      dtype=np.int32)
  Domain = domain.Domain(2,O,L,n)
  print(Domain)

  # rotation
  domain_rotation = np.deg2rad(15.0)
  Rotation = rotation.Rotation(2,domain_rotation)
  print(Rotation)

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

  GWZ = Gaussian(Domain,Rotation,ng,A,a,b,c,x0,z0)
  print(GWZ)
  
  GWZ.evaluate_gaussians()
  GWZ.plot_gaussians()
  plt.show()
  return

if __name__ == "__main__":
  test()