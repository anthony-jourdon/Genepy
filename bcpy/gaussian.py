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
    
    self.gaussian_num = None
    self.gaussian_sym = None

    domain.Domain.__init__(self,Domain.dim,Domain.O,Domain.L,Domain.n)
    rotation.Rotation.__init__(self,Domain.dim,Rotation.theta,Rotation.axis)
    return
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}:\n'
    s += f'\tDistribution for {self.ng} gaussians\n'
    for n in range(self.ng):
      s += f'\tGaussian [{n}]\n'
      s += f'\t\tAmplitude: {self.A[n]}\n'
      s += f'\t\ta: {self.a[n]}\n'
      s += f'\t\tb: {self.b[n]}\n'
      s += f'\t\tc: {self.c[n]}\n'
      s += f'\t\tCentre: [ {self.x0[n]},{self.z0[n]} ]\n'
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
    exponent = -( a*(x-x0)*(x-x0) + 2*b*(x-x0)*(z-z0) + c*(z-z0)*(z-z0) )
    if type(x) == sp.core.symbol.Symbol: u = A * sp.exp( exponent )
    else:                                u = A * np.exp( exponent )
    return u

  def symbolic_gaussian(self,A,a,b,c,x0,z0):
    g = self.gaussian_2d(A,a,b,c,self.sym_coor[0],x0,self.sym_coor[self.dim-1],z0)
    return g
  
  def numerical_gaussian(self,A,a,b,c,x0,z0):
    g = self.gaussian_2d(A,a,b,c,self.num_coor[0],x0,self.num_coor[self.dim-1],z0)
    return g
  
  def evaluate_gaussians(self):
    self.gaussian_sym = np.zeros(shape=(self.ng), dtype=object)
    self.gaussian_num = np.zeros(shape=(self.ng,*self.n), dtype=np.float64)

    for n in range(self.ng):
      if self.dim == 2:   g_centre = np.array([[self.x0[n],self.z0[n]]],dtype=np.float64)
      elif self.dim == 3: g_centre = np.array([[self.x0[n],0.0,self.z0[n]]],dtype=np.float64)
      g_centre = self.rotate_referential(g_centre,self.O,self.L)
      self.x0[n] = g_centre[0,0]
      self.z0[n] = g_centre[0,self.dim-1]
      self.gaussian_sym[n] = self.symbolic_gaussian(self.A[n],self.a[n],self.b[n],self.c[n],self.x0[n],self.z0[n])
      self.gaussian_num[n] = self.numerical_gaussian(self.A[n],self.a[n],self.b[n],self.c[n],self.x0[n],self.z0[n])
    return 
  
  def report_symbolic_functions(self):
    #s = super().report_symbolic_functions()
    s = f"Symbolic gaussian functions:\n"
    for n in range(self.ng):
      s += f"\tGaussian [{n}]:\n"
      s += f"\t\tCentre:   [ {self.x0[n]},{self.z0[n]} ]\n"
      s += f"\t\tEquation: {self.gaussian_sym[n]}\n"
    return s
  
  def compute_field_distribution(self):
    if self.gaussian_num is None:
      self.evaluate_gaussians()
    field = np.zeros(shape=(self.n), dtype=np.float64)
    for n in range(self.ng):
      field += self.gaussian_num[n]
    field = np.reshape(field,self.nv,order='F')
    return field
  
  def plot_gaussians(self):
    _, ax = plt.subplots()
    field = np.zeros(shape=(self.n), dtype=np.float64)
    for n in range(self.ng):
      field += self.gaussian_num[n] 
    if self.dim == 2:   g = ax.contourf(*self.num_coor,field,100,cmap='magma')
    elif self.dim == 3: g = ax.contourf(self.num_coor[0][:,0,:],self.num_coor[2][:,0,:],field[:,0,:],100,cmap='magma')
    ax.axis('equal')
    ax.set_title(f'Gaussian distribution {self.dim}D Domain')
    ax.set_xlabel('x axis')
    ax.set_ylabel(f'{self.sym_coor[self.dim-1]} axis')
    plt.colorbar(g,ax=ax)
    plt.draw()
    return
  
  def sprint_option(self, model_name: str):
    s = super().sprint_option(model_name)
    prefix = "wz"
    s += f"########### Initial plastic strain for weak zone ###########\n"
    s += f"-{model_name}_{prefix}_nwz {self.ng} # number of gaussians\n"
    for n in range(self.ng):
      gaussian_expression = str(self.gaussian_sym[n]).split()
      s += f"-{model_name}_{prefix}_expression_{n} "
      for term in range(len(gaussian_expression)-1):
        s += f"{gaussian_expression[term]}"
      s += f"{gaussian_expression[-1]}\n"
    return s
  
def test():
  from bcpy import utils
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
  return

def test3d():
  from bcpy import utils
  from bcpy import writers
  # Domain
  O = np.array([ 0.0, -250e3, 0.0 ], dtype=np.float64)
  L = np.array([ 600e3, 0.0, 300e3 ], dtype=np.float64)
  n = np.array([ 64, 32, 32 ], dtype=np.int32)
  Domain = domain.Domain(3,O,L,n)
  print(Domain)
  # rotation
  domain_rotation = np.deg2rad(-15.0)
  axis = np.array([0,1,0], dtype=np.float64)
  Rotation = rotation.Rotation(3,domain_rotation,axis)
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
  z0 = np.array([domain_centre[2] - dz, 
                 domain_centre[2] + dz], dtype=np.float64)
  x0[0] = utils.x_centre_from_angle(z0[0],angle,(domain_centre[0],domain_centre[2]))
  x0[1] = utils.x_centre_from_angle(z0[1],angle,(domain_centre[0],domain_centre[2]))

  GWZ = Gaussian(Domain,Rotation,ng,A,a,b,c,x0,z0)
  print(GWZ)
  GWZ.evaluate_gaussians()
  field = GWZ.compute_field_distribution()

  w = writers.WriteVTS(Domain,vtk_fname="gaussian_field.vts",point_data={"field": field})

  w.write_vts()
  GWZ.plot_gaussians()
  

if __name__ == "__main__":
  test()
  test3d()
  plt.show()