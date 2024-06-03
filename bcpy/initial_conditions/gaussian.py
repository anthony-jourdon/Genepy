#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  bc-pre-processing
#  filename: gaussian.py
#
#  This file is part of bc-pre-processing.
#
#  bc-pre-processing is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, either 
#  version 3 of the License, or any later version.
#
#  bc-pre-processing is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with bc-pre-processing. 
#  If not, see <https://www.gnu.org/licenses/>.
#====================================================================================================

import numpy as np
import sympy as sp
from bcpy.initial_conditions import domain
from bcpy import rotation

import matplotlib.pyplot as plt

class Gaussian(domain.Domain,rotation.Rotation):
  """
  .. py:class:: Gaussian(Domain,Rotation,ng,A,a,b,c,x0,z0)

    Class to build a 2D gaussian distribution defined by:
    
    .. math:: 
      u = A \\exp\\left( -\\left( a(x-x_0)^2 + 2b(x-x_0)(z-z_0) + c(z-z_0)^2 \\right) \\right)

    The class inherits from :py:class:`Domain <bcpy.Domain>` 
    and :py:class:`Rotation <bcpy.Rotation>`.

    

    :param Domain Domain: instance of the Domain class
    :param Rotation Rotation: instance of the Rotation class
    :param int ng: number of gaussians
    :param np.ndarray A: amplitude of the gaussian, shape: ``(ng,)``
    :param np.ndarray a: gaussians coefficient, shape: ``(ng,)``
    :param np.ndarray b: gaussian coefficient, shape: ``(ng,)``
    :param np.ndarray c: gaussian coefficient, shape: ``(ng,)``
    :param np.ndarray x0: x coordinate of the gaussian centre, shape: ``(ng,)``
    :param np.ndarray z0: z coordinate of the gaussian centre, shape: ``(ng,)``

    Example
    -------
    Assuming that instances of :class:`Domain` and :class:`Rotation` classes are already created
    and that 2 gaussians are required:

    .. code-block:: python

      import numpy as np
      import bcpy as bp

      ng = np.int32(2) # number of gaussians
      A  = np.array([..., ...],dtype=np.float64)
      # shape
      a  = np.array([..., ...], dtype=np.float64)
      b  = np.array([..., ...], dtype=np.float64)
      c  = np.array([..., ...], dtype=np.float64)
      x0 = np.array([..., ...], dtype=np.float64)
      z0 = np.array([..., ...], dtype=np.float64)
      # Create instance of the Gaussian class
      g  = bp.Gaussian(Domain,Rotation,ng,A,a,b,c,x0,z0)
    
    Attributes
    ----------

    .. py:attribute:: ng
      :type: int
      :canonical: bcpy.initial_conditions.gaussian.Gaussian.ng

      Number of gaussians

    .. py:attribute:: A
      :type: np.ndarray
      :canonical: bcpy.initial_conditions.gaussian.Gaussian.A

      Amplitude of the gaussian, shape: ``(ng,)``

    .. py:attribute:: a
      :type: np.ndarray
      :canonical: bcpy.initial_conditions.gaussian.Gaussian.a

      Gaussian coefficient, shape: ``(ng,)``

    .. py:attribute:: b
      :type: np.ndarray
      :canonical: bcpy.initial_conditions.gaussian.Gaussian.b

      Gaussian coefficient, shape: ``(ng,)``

    .. py:attribute:: c
      :type: np.ndarray
      :canonical: bcpy.initial_conditions.gaussian.Gaussian.c

      Gaussian coefficient, shape: ``(ng,)``

    .. py:attribute:: x0
      :type: np.ndarray
      :canonical: bcpy.initial_conditions.gaussian.Gaussian.x0

      x coordinate of the gaussian centre, shape: ``(ng,)``

    .. py:attribute:: z0
      :type: np.ndarray
      :canonical: bcpy.initial_conditions.gaussian.Gaussian.z0

      z coordinate of the gaussian centre, shape: ``(ng,)``

    .. py:attribute:: gaussian_num
      :type: np.ndarray
      :canonical: bcpy.initial_conditions.gaussian.Gaussian.gaussian_num

      Numerical values of the gaussians, shape: ``(ng,n[0],n[1])`` or ``(ng,n[0],n[1],n[2])``

    .. py:attribute:: gaussian_sym
      :type: np.ndarray
      :canonical: bcpy.initial_conditions.gaussian.Gaussian.gaussian_sym

      Symbolic values of the gaussians, shape: ``(ng,)``

    Methods
    -------
  """
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

  def gaussian_2d(self,A:float,a:float,b:float,c:float,x,x0:float,z,z0:float):
    """
    gaussian_2d(A,a,b,c,x,x0,z,z0)
    Computes a 2D gaussian distribution such that:

    .. math:: 
      u = A \\exp\\left( -\\left( a(x-x_0)^2 + 2b(x-x_0)(z-z_0) + c(z-z_0)^2 \\right) \\right)

    :param float A: amplitude of the gaussian
    :param float a,b,c: gaussian coefficients
    :param float x0: x coordinate of the centre of the gaussian
    :param float z0: z coordinate of the centre of the gaussian
    :param x,z: coordinates can be symbolic or numerical

    :return: **u**: the field with the 2D gaussian distribution
    """
    exponent = -( a*(x-x0)*(x-x0) + 2*b*(x-x0)*(z-z0) + c*(z-z0)*(z-z0) )
    if type(x) == sp.core.symbol.Symbol: u = A * sp.exp( exponent )
    else:                                u = A * np.exp( exponent )
    return u

  def symbolic_gaussian(self,A,a,b,c,x0,z0):
    """
    symbolic_gaussian(A,a,b,c,x0,z0)
    Computes the symbolic expression of a 2D gaussian distribution.
    Calls :meth:`gaussian_2d` with symbolic coordinates.

    :param float A: amplitude of the gaussian
    :param float a,b,c: gaussian coefficients
    :param float x0: x coordinate of the centre of the gaussian
    :param float z0: z coordinate of the centre of the gaussian

    :return: **u**: the symbolic expression of the field with the 2D gaussian distribution
    """
    g = self.gaussian_2d(A,a,b,c,self.sym_coor[0],x0,self.sym_coor[self.dim-1],z0)
    return g
  
  def numerical_gaussian(self,A,a,b,c,x0,z0):
    """
    numerical_gaussian(A,a,b,c,x0,z0)
    Computes the numerical values of a 2D gaussian distribution.
    Calls :meth:`gaussian_2d` with numerical coordinates.

    :param float A: amplitude of the gaussian
    :param float a,b,c: gaussian coefficients
    :param float x0: x coordinate of the centre of the gaussian
    :param float z0: z coordinate of the centre of the gaussian

    :return: **u**: the numerical values of the field with the 2D gaussian distribution
    """
    g = self.gaussian_2d(A,a,b,c,self.num_coor[0],x0,self.num_coor[self.dim-1],z0)
    return g
  
  def evaluate_gaussians(self):
    """
    evaluate_gaussians(self)
    Evaluate the symbolic and numerical values of the gaussians.
    Calls :meth:`symbolic_gaussian` and :meth:`numerical_gaussian`.
    Attach the results to the attributes 
    :attr:`gaussian_sym <bcpy.initial_conditions.gaussian.Gaussian.gaussian_sym>` 
    and :attr:`gaussian_num <bcpy.initial_conditions.gaussian.Gaussian.gaussian_num>`.

    :return: None
    """
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
    """
    report_symbolic_functions(self)
    Return a human readable string representation of
    the symbolic gaussian functions and 
    the coordinates of their centre.

    :return: **s**: string that can be printed
    """
    s = f"Symbolic gaussian functions:\n"
    for n in range(self.ng):
      s += f"\tGaussian [{n}]:\n"
      s += f"\t\tCentre:   [ {self.x0[n]},{self.z0[n]} ]\n"
      s += f"\t\tEquation: {self.gaussian_sym[n]}\n"
    return s
  
  def compute_field_distribution(self):
    """
    compute_field_distribution(self)
    Compute the gaussian distribution of a field given the 
    parameters attached to the instance of the class :class:`Gaussian`.

    :return: **field**: the field with the gaussian distribution of the shape ``(nv,)``
             with :attr:`nv <bcpy.initial_conditions.domain.Domain.nv>` the total number of nodes
             in the domain.
    :ret type: np.ndarray
    """
    if self.gaussian_num is None:
      self.evaluate_gaussians()
    field = np.zeros(shape=(self.n), dtype=np.float64)
    for n in range(self.ng):
      field += self.gaussian_num[n]
    field = np.reshape(field,self.nv,order='F')
    return field
  
  def plot_gaussians(self):
    """
    plot_gaussians(self)
    Plot a 2D view of the gaussian distribution in the domain using `matplotlib <https://matplotlib.org/>`_.

    :return: None
    """
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
    """
    sprint_option(self,model_name:str)
    Return a string formatted for `pTatin3d`_ input file using `PETSc`_ options format.

    :param str model_name: name of the model

    :return: string with the options
    """
    prefix = "wz"
    s  = f"########### Initial plastic strain for weak zone ###########\n"
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
