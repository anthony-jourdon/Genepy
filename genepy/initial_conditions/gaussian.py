#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: gaussian.py
#
#  This file is part of Genepy.
#
#  Genepy is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, either 
#  version 3 of the License, or any later version.
#
#  Genepy is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with Genepy. 
#  If not, see <https://www.gnu.org/licenses/>.
#====================================================================================================

import numpy as np
import sympy as sp
from genepy.initial_conditions import domain
from genepy import rotation
from collections.abc import Iterable

import matplotlib.pyplot as plt

class Gaussian(domain.Domain):
  """
  .. py:class:: Gaussian(Domain)

    Parent class of all the Gaussian distributions.

    :Attributes:

    .. py:attribute:: gaussian_sym
      :type: sympy.Expr
      :canonical: genepy.initial_conditions.gaussian.Gaussian.gaussian_sym

      Symbolic expression of the gaussian distribution

    .. py:attribute:: gaussian_num
      :type: np.ndarray
      :canonical: genepy.initial_conditions.gaussian.Gaussian.gaussian_num

      Numerical values of the gaussian distribution
  """
  def __init__(self,Domain:domain.Domain) -> None:
    domain.Domain.__init__(self,Domain.dim,Domain.O_num,Domain.L_num,Domain.n)

    self.gaussian_sym = None
    self.gaussian_num = None
    return
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'Symbolic expression: {self.gaussian_sym}\n'
    return s
    

class GaussianConstructor(Gaussian):
  """
 .. py:class:: GaussianConstructor(Domain,A,a,expression)
  
    Child class of :py:class:`Gaussian <genepy.initial_conditions.gaussian.Gaussian>` 
    to build a function describing a gaussian distribution defined by:

    .. math::
      g(\\mathbf a, \\mathbf x) = A \\exp\\left( - \\mathbf a \\cdot \\mathbf x \\cdot \\mathbf x \\right)

    where :math:`\\mathbf a` is an array of coefficients such that :math:`\\mathbf a = [a_1,a_2,...,a_n]`, 
    :math:`\\mathbf x` is an array such that :math:`\\mathbf x = [x_1,x_2,...,x_n]` where
    :math:`x_i, \\, i = 1,2,...,n` can be single variables or mathematical functions and 
    :math:`A` is the amplitude of the gaussian.

    This class is low level and should be called by the user to create its own gaussian distribution if not 
    already available with the child classes.
  
    :param Domain Domain: instance of the Domain class
    :param float A: amplitude of the gaussian
    :param float | list | numpy.ndarray a: gaussian extent coefficient(s). Can be a single float, a list or a numpy array.
    :param sympy.Expr | list | numpy.ndarray expression: expression(s) of the gaussian. Can be a single expression, a list or a numpy array.
  """
  def __init__(
      self,
      Domain:domain.Domain,
      A:float,
      a:float|list[float]|np.ndarray,
      expression:sp.Expr|list[sp.Expr]|np.ndarray,
    ) -> None:
    Gaussian.__init__(self,Domain)

    self.A = A
    self.a = a
    self.x = expression
    self.gaussian()
    return
  
  def __str__(self) -> str:
    s = Gaussian.__str__(self)
    s += f'Amplitude:      {self.A}\n'
    s += f'Coefficient(s): {self.a}\n'
    s += f'Exponent(s):    {self.x}\n'
    return s

  def _gaussian(self,A,a,x):
    """
    Compute the gaussian distribution 
    
    .. math::

      g(a,x) = A \\exp\\left( -a x^2 \\right)

    for a single dimension.
    """
    return A * sp.exp( -a*x**2 )
  
  def gaussian(self):
    """
    Compute the gaussian distribution for the given expression(s).
    Make use of the rule :math:`\\exp(a+b) = \\exp(a)\\exp(b)` to build the final multidimensional expression.
    The method is called by the class constructor and does not require to be called by the user.
    """
    # check the dimensionnality of the expression
    if isinstance(self.a,Iterable):
      nd = len(self.a) # number of dimensions
      g = self._gaussian(self.A,self.a[0],self.x[0]) # first dimension
      for n in range(1,nd):
        g *= self._gaussian(self.A,self.a[n],self.x[n]) # next dimensions
    else:
      g = self._gaussian(self.A,self.a,self.x)
    self.gaussian_sym = sp.powsimp(g) # simplify the expression
    return
  
  def evaluate_gaussian(self):
    """
    Evaluate numerically the expression of the gaussian distribution.
    Uses the method :py:meth:`sympy.lambdify` to convert the symbolic expression to a python callable.
    """
    if self.gaussian_sym is None:
      self.gaussian()
    g = sp.lambdify(self.sym_coor,self.gaussian_sym,"numpy")
    self.gaussian_num = g(*self.num_coor)
    return
    
class Gaussian2D(GaussianConstructor):
  """
  .. py:class:: Gaussian2D(Domain,A,a,b,x0,z0,Rotation=None)

    Child class of :py:class:`GaussianConstructor <genepy.initial_conditions.gaussian.GaussianConstructor>` 
    to build a 2D gaussian distribution defined by:

    .. math::

      g((a,b);(x,z)) = A \\exp\\left( -\\left( a(x-x_0)^2 + b(z-z_0)^2 \\right) \\right)
    
    where :math:`a` and :math:`b` are coefficients controlling the shape of the gaussian in 
    the :math:`x` and :math:`z` directions respectively, 
    :math:`x_0` and :math:`z_0` are the coordinates of the centre of the gaussian 
    and :math:`A` is the amplitude of the gaussian.

    :param Domain Domain: instance of the Domain class
    :param float A: amplitude of the gaussian
    :param float a: gaussian extent coefficient in the :math:`x` direction
    :param float b: gaussian extent coefficient in the :math:`z` direction
    :param float x0: :math:`x` coordinate of the centre of the gaussian
    :param float z0: :math:`z` coordinate of the centre of the gaussian
    :param Rotation Rotation: instance of the Rotation class (**optional**) to rotate the centre of the gaussian
    
    :Example:

    .. code:: python

      import numpy as np
      import genepy as gp

      # Domain
      O = np.array([ 0.0, -250e3, 0.0 ], dtype=np.float64)
      L = np.array([ 600e3, 0.0, 300e3 ], dtype=np.float64)
      n = np.array([ 64, 32, 64 ], dtype=np.int32)
      Domain = gp.Domain(3,O,L,n)

      a = 0.5 * 6.0e-5**2
      b = 0.5 * 6.0e-5**2
      x0 = 0.5*Domain.L_num[0]
      z0 = 0.5*Domain.L_num[2]

      Gaussian = gp.Gaussian2D(Domain,1.0,a,b,x0,z0)
      print(Gaussian) # prints the expression and the parameters of the gaussian
    
    If one wants to numerically evaluate the gaussian distribution use:

    .. code:: python

      Gaussian.evaluate_gaussian()

  """
  def __init__(
      self,
      Domain:domain.Domain,
      A:float, # amplitude of the gaussian
      a:float, # gaussian extent coefficient in the x direction
      b:float, # gaussian extent coefficient in the z direction
      x0:float, # x coordinate of the centre of the gaussian
      z0:float, # z coordinate of the centre of the gaussian
      Rotation:rotation.Rotation|None=None
    ) -> None:

    if Rotation is not None:
      if Domain.dim == 2:   g_centre = np.array([ [x0, z0] ],    dtype=np.float64)
      elif Domain.dim == 3: g_centre = np.array([ [x0, 0, z0] ], dtype=np.float64)
      g_centre = Rotation.rotate_referential(g_centre,Domain.O_num,Domain.L_num)
      x0 = g_centre[ 0, 0 ]
      z0 = g_centre[ 0, Domain.dim-1 ]

    expression = [
      Domain.sym_coor[0] - x0,
      Domain.sym_coor[2] - z0
    ]
    coefficients = [a,b]
    GaussianConstructor.__init__(self,Domain,A,coefficients,expression)
    return

class Gaussian3D(GaussianConstructor):
  """
  .. py:class:: Gaussian3D(Domain,A,a,b,c,x0,y0,z0,Rotation=None)

    Child class of :py:class:`GaussianConstructor <genepy.initial_conditions.gaussian.GaussianConstructor>` 
    to build a 3D gaussian distribution defined by:

    .. math::

      g((a,b,c);(x,y,z)) = A \\exp\\left( -\\left( a(x-x_0)^2 + b(y-y_0)^2 + c(z-z_0)^2 \\right) \\right)
    
    where :math:`a`, :math:`b` and :math:`c` are coefficients controlling the shape of the gaussian in 
    the :math:`x`, :math:`y` and :math:`z` directions respectively, 
    :math:`x_0`, :math:`y_0` and :math:`z_0` are the coordinates of the centre of the gaussian 
    and :math:`A` is the amplitude of the gaussian.

    :param Domain Domain: instance of the Domain class
    :param float A: amplitude of the gaussian
    :param float a: gaussian extent coefficient in the :math:`x` direction
    :param float b: gaussian extent coefficient in the :math:`y` direction
    :param float c: gaussian extent coefficient in the :math:`z` direction
    :param float x0: :math:`x` coordinate of the centre of the gaussian
    :param float y0: :math:`y` coordinate of the centre of the gaussian
    :param float z0: :math:`z` coordinate of the centre of the gaussian
    :param Rotation Rotation: instance of the Rotation class (**optional**) to rotate the centre of the gaussian
    
    :Example:

    .. code:: python

      import numpy as np
      import genepy as gp

      # Domain
      O = np.array([ 0.0, -250e3, 0.0 ], dtype=np.float64)
      L = np.array([ 600e3, 0.0, 300e3 ], dtype=np.float64)
      n = np.array([ 64, 32, 64 ], dtype=np.int32)
      Domain = gp.Domain(3,O,L,n)

      a = 0.5 * 6.0e-5**2
      b = 0.5 * 6.0e-5**2
      c = 0.5 * 6.0e-5**2
      x0 = 0.5*Domain.L_num[0]
      y0 = 0.5*Domain.L_num[1]
      z0 = 0.5*Domain.L_num[2]

      Gaussian = gp.Gaussian3D(Domain,1.0,a,b,c,x0,y0,z0)
      print(Gaussian) # prints the expression and the parameters of the gaussian
    
    If one wants to numerically evaluate the gaussian distribution use:

    .. code:: python

      Gaussian.evaluate_gaussian()

  """
  def __init__(
    self, 
    Domain: domain.Domain, 
    A: float, # amplitude of the gaussian
    a: float, # gaussian extent coefficient in the x direction
    b: float, # gaussian extent coefficient in the y direction
    c: float, # gaussian extent coefficient in the z direction
    x0: float, # x coordinate of the centre of the gaussian
    y0: float, # y coordinate of the centre of the gaussian
    z0: float, # z coordinate of the centre of the gaussian
    Rotation: rotation.Rotation | None = None
  ) -> None:

    if Rotation is not None:
      g_centre = np.array([ [x0, y0, z0] ], dtype=np.float64)
      g_centre = Rotation.rotate_referential(g_centre,Domain.O_num,Domain.L_num)
      x0 = g_centre[ 0, 0 ]
      y0 = g_centre[ 0, 1 ]
      z0 = g_centre[ 0, 2 ]

    expression = [
      Domain.sym_coor[0] - x0,
      Domain.sym_coor[1] - y0,
      Domain.sym_coor[2] - z0
    ]
    coefficients = [a,b,c]
    GaussianConstructor.__init__(self,Domain,A,coefficients,expression)
    return

class GaussianPlane(GaussianConstructor):
  """
  .. py:class:: GaussianPlane(Domain,A,a,plane_coeff)

    Child class of :py:class:`GaussianConstructor <genepy.initial_conditions.gaussian.GaussianConstructor>`
    to build a 3D gaussian distribution around a plane defined by:

    .. math::

      g( a, D(\\mathbf x) ) = A \\exp\\left( - a D(\\mathbf x)^2 \\right)

    where :math:`a` is the coefficient controlling the shape of the gaussian around the plane and 

    .. math::

      D(\\mathbf x) = \\frac{n_0 x + n_1 y + n_2 z + d}{\\sqrt{n_0^2 + n_1^2 + n_2^2}} 
      = \\frac{\\mathbf n \\cdot \\mathbf x + d}{|| \\mathbf n ||} 

    is the equation describing the distance between a point of coordinates :math:`\\mathbf x` 
    and the plane defined by the normal vector :math:`\\mathbf n = [n_0, n_1, n_2]` and the parameter 
    :math:`d`.

    :param Domain Domain: instance of the Domain class
    :param float A: amplitude of the gaussian
    :param float a: gaussian extent coefficient
    :param np.ndarray plane_coeff: plane coefficients: ``numpy.array([n_0, n_1, n_2, d])``
                                   where :math:`n_0 x + n_1 y + n_2 z + d = 0`
    
    :Example:

    .. code:: python

      import numpy as np
      import genepy as gp

      # Domain
      O = np.array([ 0.0, -250e3, 0.0 ], dtype=np.float64)
      L = np.array([ 600e3, 0.0, 300e3 ], dtype=np.float64)
      n = np.array([ 64, 32, 64 ], dtype=np.int32)
      Domain = gp.Domain(3,O,L,n)

      # gaussian extent coefficient
      a = 0.5 * 6.0e-5**2

      # points defining the plane
      pt_A = np.array([200.0, 0.0, 0.0], dtype=np.float64) * 1e3
      pt_B = np.array([200.0, 0.0, 300.0], dtype=np.float64) * 1e3
      pt_C = np.array([600.0, -100.0, 0.0], dtype=np.float64) * 1e3
      
      # normal vector to the plane
      normal = np.cross(pt_B - pt_A, pt_C - pt_A)
      # d parameter
      d = -np.dot(normal,pt_A)

      plane_coeff = np.array([normal[0], normal[1], normal[2], d], dtype=np.float64)

      Gaussian = gp.GaussianPlane(Domain,1.0,a,plane_coeff)

  """
  def __init__(
      self,
      Domain:domain.Domain,
      A:float, # amplitude of the gaussian
      a:float, # gaussian extent coefficient
      coeff:np.ndarray, # [a,b,c,d] -> plane coefficients: a*x + b*y + c*z + d
    ) -> None:

    if isinstance(coeff,Iterable):
      if len(coeff) != 4:
        raise RuntimeError(f"Plane coefficients must have 4 elements, found {len(coeff)}")
    else:
      raise RuntimeError(f"Plane coefficients must be an iterable, found {type(coeff)}")

    expression = (
      coeff[0]*Domain.sym_coor[0] + # a*x
      coeff[1]*Domain.sym_coor[1] + # b*y
      coeff[2]*Domain.sym_coor[2] + # c*z
      coeff[3]                      # d
    ) / sp.sqrt(coeff[0]**2 + coeff[1]**2 + coeff[2]**2) # norm of normal vector
    GaussianConstructor.__init__(self,Domain,A,a,expression)
    return
  
class GaussianCircle(GaussianConstructor):
  """
  .. py:class:: GaussianCircle(Domain,A,a,xc,zc,r,Rotation=None)

    Child class of :py:class:`GaussianConstructor <genepy.initial_conditions.gaussian.GaussianConstructor>`
    to build a 2D gaussian distribution around a circle defined by:

    .. math::

      g( a, D(\\mathbf x) ) = A \\exp\\left( - a D(\\mathbf x)^2 \\right)

    where :math:`a` is the coefficient controlling the shape of the gaussian around the circle and

    .. math::

      D(\\mathbf x) = \\sqrt{(x-x_c)^2 + (z-z_c)^2} - r

    is the equation describing the distance between a point of coordinates :math:`\\mathbf x`
    and the circle defined by the centre :math:`(x_c,z_c)` and the radius :math:`r`.

    :param Domain Domain: instance of the Domain class
    :param float A: amplitude of the gaussian
    :param float a: gaussian extent coefficient
    :param float xc: :math:`x` coordinate of the centre of the circle
    :param float zc: :math:`z` coordinate of the centre of the circle
    :param float r: radius of the circle
    :param Rotation Rotation: instance of the Rotation class (**optional**) to rotate the centre of the circle

    :Example:

    .. code:: python

      import numpy as np
      import genepy as gp

      # Domain
      O = np.array([ 0.0, -250e3, 0.0 ], dtype=np.float64)
      L = np.array([ 600e3, 0.0, 300e3 ], dtype=np.float64)
      n = np.array([ 64, 32, 64 ], dtype=np.int32)
      Domain = gp.Domain(3,O,L,n)

      # gaussian extent coefficient
      coeff = 0.5*6.0e-5**2
      # centre of the circle
      xc = 0.5*Domain.L_num[0]
      zc = 0.5*Domain.L_num[2]
      # radius of the circle
      r  = 0.25*Domain.L_num[2]

      Gaussian = gp.GaussianCircle(Domain,1.0,coeff,xc,zc,r)
      
  """
  def __init__(
      self,
      Domain:domain.Domain,
      A:float, # amplitude of the gaussian
      a:float, # gaussian extent coefficient
      xc:float, # x coordinate of the centre of the circle
      zc:float, # z coordinate of the centre of the circle
      r:float,  # radius of the circle
      Rotation:rotation.Rotation|None=None
  ) -> None:
    
    if Rotation is not None:
      if Domain.dim == 2:   g_centre = np.array([ [xc, zc] ],    dtype=np.float64)
      elif Domain.dim == 3: g_centre = np.array([ [xc, 0, zc] ], dtype=np.float64)
      g_centre = Rotation.rotate_referential(g_centre,Domain.O_num,Domain.L_num)
      xc = g_centre[ 0, 0 ]
      zc = g_centre[ 0, Domain.dim-1 ]

    expression = sp.sqrt((Domain.sym_coor[0]-xc)**2 + (Domain.sym_coor[2]-zc)**2) - r
    GaussianConstructor.__init__(self,Domain,A,a,expression)
    return

class GaussiansOptions:
  """
  .. py:class:: GaussiansOptions(gaussians,blocksize=10)

    Class to generate the options for the input file of `pTatin3d`_.
    The class takes a list of gaussians and can sum them in blocks of size `blocksize` 
    to reduce the number of options.
    The class instance can be used to generate the options for the gaussian distributions.
    If the gaussians are used to set an initial plastic strain, it should be passed to the class
    :py:class:`InitialPlasticStrain <genepy.initial_conditions.plastic_strain.InitialPlasticStrain>`.

    :param list gaussians: list of Gaussian instances
    :param int blocksize: (**Optional**) size of the block to sum the gaussians.

    :Example:

    .. code:: python

      import numpy as np
      import genepy as gp

      # Domain
      O = np.array([ 0.0, -250e3, 0.0 ], dtype=np.float64)
      L = np.array([ 600e3, 0.0, 300e3 ], dtype=np.float64)
      n = np.array([ 64, 32, 64 ], dtype=np.int32)
      Domain = gp.Domain(3,O,L,n)

      # gaussian distribution
      coeff = 0.5*6.0e-5**2
      x0 = [0.25 * Domain.L_num[0], 0.75 * Domain.L_num[0]]
      z0 = [0.25 * Domain.L_num[2], 0.75 * Domain.L_num[2]]
      gaussians = []
      for i in range(2):
        gaussians.append(gp.Gaussian2D(Domain,1.0,coeff,coeff,x0[i],z0[i]))
      
      Gopt = gp.GaussiansOptions(gaussians,blocksize=5)

    Pass the instance ``Gopts`` to the class 
    :py:class:`InitialPlasticStrain <genepy.initial_conditions.plastic_strain.InitialPlasticStrain>`:

    .. code:: python

      IPS = gp.InitialPlasticStrain(Gopt)
      
  """
  def __init__(self,gaussians:list[Gaussian],blocksize:int=10) -> None:
    self.gaussians = gaussians
    self.ng        = len(gaussians)
    self.gaussian_sym = []
    if len(gaussians) <= blocksize:
      for g in gaussians:
        self.gaussian_sym.append(g.gaussian_sym)
    else:
      self.sum_gaussians(blocksize)
    return
  
  def sum_gaussians(self,blocksize):
    """
    sum_gaussians(self,blocksize)
    Sum the gaussians in blocks of size `blocksize`.
    Can be utilized to reduce the amount of options in the input file for `pTatin3d`_.

    :param int blocksize: size of the block to sum the gaussians.
    """
    nblocks = int(np.ceil(self.ng/blocksize))
    for n in range(nblocks):
      start = n*blocksize
      end   = min((n+1)*blocksize,self.ng)
      for i in range(start,end):
        if i == start: field_sym  = self.gaussians[i].gaussian_sym
        else:          field_sym += self.gaussians[i].gaussian_sym
      self.gaussian_sym.append(field_sym)
    self.ng = nblocks
    return 

  def sprint_option(self, model_name:str, prefix:str) -> str:
    if prefix == "wz":
      s = f"-{model_name}_{prefix}_nwz {self.ng} # number of gaussians\n"
    elif prefix == "heat_source":
      s = f"-{model_name}_{prefix}_nhs {self.ng} # number of gaussians\n"    
    for n in range(self.ng):
      gaussian_expression = str(sp.ccode(self.gaussian_sym[n])).replace(" ","")
      s += f"-{model_name}_{prefix}_expression_{n} {gaussian_expression}\n"
    return s

def test():
  from genepy import utils
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

  GWZ = Gaussian(Domain,ng,A,a,b,c,x0,z0,Rotation)
  print(GWZ)
  
  GWZ.evaluate_gaussians()
  GWZ.plot_gaussians()
  return

def test3d():
  from genepy import utils
  from genepy import writers
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

  GWZ = Gaussian(Domain,ng,A,a,b,c,x0,z0,Rotation)
  print(GWZ)
  GWZ.evaluate_gaussians()
  field = GWZ.compute_field_distribution()

  w = writers.WriteVTS(Domain,vtk_fname="gaussian_field.vts",point_data={"field": field})

  w.write_vts()
  GWZ.plot_gaussians()
  
def test_function_gaussian():
  from collections.abc import Iterable
  # Domain
  O = np.array([ 0.0, -250e3, 0.0 ], dtype=np.float64)
  L = np.array([ 600e3, 0.0, 300e3 ], dtype=np.float64)
  n = np.array([ 64, 32, 32 ], dtype=np.int32)
  Domain = domain.Domain(3,O,L,n)
  
  # standard 2D gaussian
  f1 = Domain.sym_coor[0] - 0.5*Domain.L_num[0]
  f2 = Domain.sym_coor[2] - 0.5*Domain.L_num[2]
  print(isinstance(f1,Iterable))
  print(isinstance(f1,sp.Expr))

  a1 = 0.5 * 6.0e-5**2
  a2 = 0.5 * 3.0e-5**2
  g1 = sp.exp(-a1*f1**2)
  g2 = sp.exp(-a2*f2**2)
  G = sp.powsimp(g1*g2)
  GG = sp.lambdify((Domain.sym_coor[0],Domain.sym_coor[2]),G,"numpy")
  H = GG(Domain.num_coor[0],Domain.num_coor[2])

  # funny gaussian distribution around a circle
  C = np.array([0.5*Domain.L_num[0], 0.5*Domain.L_num[2]],dtype=np.float64)
  r = 0.25*Domain.L_num[2]
  f = sp.sqrt((Domain.sym_coor[0]-C[0])**2 + (Domain.sym_coor[2]-C[1])**2) - r
  
  g = sp.exp(-a1*f**2)
  G = sp.powsimp(g)
  GG = sp.lambdify((Domain.sym_coor[0],Domain.sym_coor[2]),G,"numpy")
  F = GG(Domain.num_coor[0],Domain.num_coor[2])

  _, ax = plt.subplots(ncols=1,nrows=2,tight_layout=True)
  c = ax[0].contourf(Domain.num_coor[0][:,0,:],Domain.num_coor[2][:,0,:],H[:,0,:],100,cmap='magma')
  ax[0].axis('equal')
  ax[0].set_title(f'2D Gaussian distribution')
  ax[0].set_xlabel('x axis')
  ax[0].set_ylabel('z axis')

  c = ax[1].contourf(Domain.num_coor[0][:,0,:],Domain.num_coor[2][:,0,:],F[:,0,:],100,cmap='magma')
  ax[1].axis('equal')
  ax[1].set_title(f'Gaussian distribution around a circle')
  ax[1].set_xlabel('x axis')
  ax[1].set_ylabel('z axis')
  
  plt.draw()
  plt.show()
  return

def test_single_gaussian_class():
  # Domain
  O = np.array([ 0.0, -250e3, 0.0 ], dtype=np.float64)
  L = np.array([ 600e3, 0.0, 300e3 ], dtype=np.float64)
  n = np.array([ 64, 32, 32 ], dtype=np.int32)
  Domain = domain.Domain(3,O,L,n)

  # standard 2D gaussian
  """
  expression = [
    Domain.sym_coor[0] - 0.5*Domain.L_num[0],
    Domain.sym_coor[2] - 0.5*Domain.L_num[2]
  ]
  a = [
    0.5 * 6.0e-5**2,
    0.5 * 3.0e-5**2
  ]
  """
  C = np.array([0.5*Domain.L_num[0], 0.5*Domain.L_num[2]],dtype=np.float64)
  r = 0.25*Domain.L_num[2]
  expression = sp.sqrt((Domain.sym_coor[0]-C[0])**2 + (Domain.sym_coor[2]-C[1])**2) - r
  a = 0.5 * 6.0e-5**2
  G = GaussianConstructor(Domain,1.0,a,expression)
  print(G.gaussian_sym)
  G.evaluate_gaussian()
  _,ax = plt.subplots()
  c = ax.contourf(Domain.num_coor[0][:,0,:],Domain.num_coor[2][:,0,:],G.gaussian_num[:,0,:],100,cmap='magma')
  ax.axis('equal')
  ax.set_title(f'2D Gaussian distribution')
  ax.set_xlabel('x axis')
  ax.set_ylabel('z axis')
  plt.draw()
  plt.show()
  return

if __name__ == "__main__":
  #test()
  #test3d()
  #plt.show()
  #test_function_gaussian()
  test_single_gaussian_class()
