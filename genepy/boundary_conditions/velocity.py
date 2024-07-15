#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: velocity.py
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

from genepy.initial_conditions import domain
from genepy import rotation
from genepy import writers
from genepy.utils import newton_raphson
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

class Velocity(domain.Domain,rotation.Rotation):
  """
  .. py:class:: Velocity(Domain, Rotation=None)

    Class to evaluate symbolic and numeric velocity function and its gradient.
    This class is the parent class of the different velocity classes and contains methods shared by all.
    Case specific methods are implemented in the corresponding child classes.
    The class inherits from :py:class:`Domain <genepy.Domain>` and :py:class:`Rotation <genepy.Rotation>` classes.

    Methods
    -------
  """
  def __init__(self,Domain:domain.Domain,Rotation:rotation.Rotation=None) -> None:
    if isinstance(Domain, domain.DomainALE):
      domain.DomainALE.__init__(self, Domain.dim, Domain.O_num, Domain.L_num, Domain.n, coor=Domain.num_coor)
    else:
      domain.Domain.__init__(self, Domain.dim, Domain.O, Domain.L, Domain.n, coor=Domain.num_coor)
    # set rotation angle to zero if Rotation class is not provided
    if Rotation is None: rotation.Rotation.__init__(self,Domain.dim,0.0,np.array([0,1,0]))
    else:                rotation.Rotation.__init__(self,Domain.dim,Rotation.theta,Rotation.axis)
    self.u                = None
    self.grad_u           = None
    self.uO               = None
    self.uL               = None
    self.u_dir_horizontal = None
    self.u_dir            = None
    return

  def report_symbolic_functions(self) -> str:
    """
    report_symbolic_functions(self)
    Returns a string with the symbolic velocity function, 
    its gradient and the boundary velocity orientation.
    Can be used with print() or written to a file.

    :return: string with the symbolic velocity function, its gradient and the boundary velocity orientation.
    :rtype: str

    """
    s = f"Symbolic velocity function:\n"
    for i in range(self.dim):
      s += f"\tu{self.sym_coor[i]}{self.sym_coor} = {self.u[i]}\n"
    if self.grad_u is not None:
      s += f"Gradient of the velocity function:\n"
      for i in range(self.dim):
        for j in range(self.dim):
          s += f"\tdu{self.sym_coor[i]}/d{self.sym_coor[j]} = "+str(self.grad_u[i,j])+"\n"
    if self.u_dir_horizontal is not None and self.u_dir is not None:
      s += f"Boundary velocity orientation (unit vector):\n"
      s += "\tHorizontal vector:\n"
      for i in range(self.dim):
        s += f"\t\tu_dir_{self.sym_coor[i]} = {self.u_dir_horizontal[i]}\n"
      s += "\tFull vector:\n"
      for i in range(self.dim):
        s += f"\t\tu_dir_{self.sym_coor[i]} = {self.u_dir[i]}\n"
    return s
  
  def velocity_function(self):
    s = "velocity_function(self)\n"
    s += "This method is implemented in the child classes.\n"
    s += "Use one of the child classes to evaluate the velocity field."
    raise NotImplementedError(s)

  def evaluate_gradient(self,u):
    """
    evaluate_gradient(self,u)
    Evaluates the gradient of the velocity function in symbolic form 
    returning a matrix of the shape ``(dim,dim)`` such that:
      
    .. math:: 
      \\left( \\nabla \\mathbf u \\right)_{ij} = \\frac{\\partial u_i}{\\partial x_j}

    
    :param u: vector valued function of the velocity field

    :return: **grad_u**: matrix of the gradient of the velocity field shape ``(dim,dim)``
    """
    grad_u = np.zeros(shape=(self.dim,self.dim), dtype='object')
    for i in range(self.dim):
      for j in range(self.dim):
        grad_u[i,j] = u[i].diff(self.sym_coor[j])
    return grad_u
  
  def evaluate_u_dot_n(self,u):
    """
    evaluate_u_dot_n(self,u)
    Evaluates the dot product of the velocity field with the normal vector in symbolic form 
    such that

    .. math:: 
      \\mathbf u \\cdot \\mathbf n = \\sum_{i=1}^{d} u_i n_i
    
    with :math:`\\mathbf u` the vector valued velocity function, 
    :math:`\\mathbf n` the normal vector to the boundary pointing outward the domain and
    :math:`d` the number of spatial dimensions.  
    
    :param u: Vector valued function of the velocity field.
              If not provided, the symbolic velocity field is evaluated with :meth:`evaluate_velocity_symbolic`
    
    :return: **u_dot_n**: dot product of the velocity field with the normal vector
    """
    nmap = {1:'n_x',2:'n_x n_y',3:'n_x n_y n_z'}
    n = sp.symbols(nmap[self.dim])
    
    u = sp.Matrix(u)
    n = sp.Matrix([n])
    u_dot_n = u.dot(n)
    return u_dot_n

  def evaluate_int_u_dot_n_faces(self,u):
    """
    evaluate_int_u_dot_n_faces(self,u)
    Evaluates the integral of the dot product of the velocity field 
    with the normal vector over the faces of the domain such that

    .. math::
      I = \\int_S \\mathbf u \\cdot \\mathbf n \\, dS

    The integral is computed over each face of the domain and stored in a dictionary:

    .. code-block:: python
      
        int_u_dot_n = {
          'xmin': int_u_dot_n_dxmin,
          'xmax': int_u_dot_n_dxmax,
          'ymin': int_u_dot_n_dymin,
          'ymax': int_u_dot_n_dymax,
          'zmin': int_u_dot_n_dzmin,
          'zmax': int_u_dot_n_dzmax
        }
    
    :param u: Vector valued function of the velocity field.
              If not provided, the symbolic velocity field is evaluated with :meth:`evaluate_velocity_symbolic`

    :return: **int_u_dot_n**: dictionary with the integral of the dot product of the velocity field with the normal vector over the faces of the domain
    :rtype: dict
    """
    # evaluate symbolic u.n
    u_dot_n = self.evaluate_u_dot_n(u)
    directions = ["x","y","z"]
    faces      = ["min","max"]
    
    int_u_dot_n = {}
    for d in range(self.dim): # loop over each direction
      for f in range(2): # loop over faces (min,max)
        # initialize normal to zero
        normal = np.zeros(shape=(self.dim), dtype=np.float64)
        # construct face name
        face = directions[d]+faces[f]
        # set non zero component of the normal vector -1 if min face, 1 if max face 
        normal[d] = (-1)**(f+1)
        # substitute the normal vector components by numerical values
        n_num = {'n_x':normal[0],'n_y':normal[1]}
        if self.dim == 3: n_num['n_z'] = normal[2]
        # evaluate u.n
        udn = u_dot_n.subs(n_num)

        # substitute the coordinate components by values corresponding to the face
        if type(self.O[0]) is sp.Symbol: dtype = 'object'
        else:                            dtype = np.float64
        _x = np.zeros(shape=(self.dim), dtype=dtype)
        if faces[f] == "min": _x[d] = self.O[d]
        if faces[f] == "max": _x[d] = self.L[d]
        # evaluate u.n at the face
        udn = udn.subs({self.sym_coor[d]:_x[d]})

        # integrate u.n over the face: simple integral if 2d, double integral if 3d
        if d == 0:
          # I2 = int_Oy^Ly u.n dy
          integral_u_dot_n = sp.integrate(udn,(self.sym_coor[1],self.O[1],self.L[1]))
          if self.dim == 3: 
            # I3 = int_Oz^Lz I2 dz
            integral_u_dot_n = sp.integrate(integral_u_dot_n,(self.sym_coor[2],self.O[2],self.L[2]))
        if d == 1:
          integral_u_dot_n = sp.integrate(udn,(self.sym_coor[0],self.O[0],self.L[0]))
          if self.dim == 3: 
            integral_u_dot_n = sp.integrate(integral_u_dot_n,(self.sym_coor[2],self.O[2],self.L[2]))
        if d == 2:
          integral_u_dot_n = sp.integrate(udn,(self.sym_coor[0],self.O[0],self.L[0]))
          integral_u_dot_n = sp.integrate(integral_u_dot_n,(self.sym_coor[1],self.O[1],self.L[1]))
        

        int_u_dot_n[face] = integral_u_dot_n.simplify()
    return int_u_dot_n

  def velocity_boundary(self):
    if self.u is None: raise RuntimeError("The velocity field has not been evaluated yet.")
    self.uO = np.zeros(shape=(self.dim), dtype='object')
    self.uL = np.zeros(shape=(self.dim), dtype='object')
    coor = {self.sym_coor[d]:self.O_num[d] for d in range(self.dim)}
    for d in range(self.dim):
      self.uO[d] = self.u[d].subs(coor)
      self.uL[d] = self.u[d].subs(coor)
    print("Boundary velocity at O:",self.uO)
    print("Boundary velocity at L:",self.uL)
    return

  def get_velocity_orientation(self,horizontal=True,normalize=True):
    """
    get_velocity_orientation(self,horizontal=True,normalize=False)
    Returns the orientation vector of the velocity field at the boundary

    :param bool horizontal: if True, only the horizontal components are returned (default: True)
    :param bool normalize:  if True, the vector is normalized (default: True)

    :return: **uL**: orientation of the velocity field at the boundary
    """
    if self.uL is None or self.uO is None: self.velocity_boundary()
    if horizontal == True:
      uL = np.zeros(shape=(self.dim), dtype=np.float64)
      uL[0] = self.uL[0]
      uL[self.dim-1] = self.uL[self.dim-1]
    else:
      uL = np.copy(self.uL)
    R  = self.rotation_matrix()
    uL = self.rotate_vector(R,uL,ccw=True)
    if normalize == True: 
      uL = uL / np.linalg.norm(uL)
    return uL

class VelocityLinear(Velocity):
  """
  .. py:class:: VelocityLinear(Domain, u_norm, variation_dir, velocity_type, u_angle=0.0, Rotation=None)

    Class to evaluate symbolic and numeric linear velocity function and its gradient in space.
    The velocity field is defined as a linear function of the coordinates. 
    It can be used for 2D and 3D domains.
    Given a horizontal velocity field, the vertical velocity is computed by integrating the dot product of the velocity field with the normal vector over the faces of the domain.
    The class inherits from :py:class:`Domain <genepy.initial_conditions.domain.Domain>` and :py:class:`Rotation <genepy.rotation.Rotation>`

    :param Domain Domain: domain in which the velocity field is evaluated
    :param float u_norm: velocity norm of the vector along boundaries
    :param str variation_dir: direction in which the velocity varies (``"x"``, ``"y"``, ``"z"``)
    :param str velocity_type: velocity field orientation, 
                              ``"extension"`` velocity is directed outward the domain, 
                              ``"compression"`` velocity is directed inward the domain
    
    **Optional**
    
    :param float u_angle: angle in radians of the velocity field with the z axis
    :param Rotation Rotation: Rotation object to rotate the referential

    Examples
    --------
    Assuming the class :class:`Domain <genepy.initial_conditions.domain.Domain>` as been instanciated as ``Domain``
    
    **Without rotation** of the referential

    .. code-block:: python

      u_norm  = 1.0              # horizontal velocity norm
      u_angle = np.deg2rad(45.0) # velocity angle in [-pi/2, pi/2]. If 0, can be ommited
      u_dir   = "z"              # direction in which velocity varies
      u_type  = "extension"      # velocity orientation
      # Create VelocityLinear instance 
      Vel = genepy.VelocityLinear(Domain,u_norm,u_dir,u_type,u_angle)
    
    **With rotation** of the referential

    .. code-block::  python

      # Rotation of the referential
      dim     = 3                                   # Number of spatial dimensions
      r_angle = np.deg2rad(-15.0)                   # Rotation angle
      axis    = np.array([0,1,0], dtype=np.float64) # Rotation axis
      # Create Rotation instance
      Rotation = gp.Rotation(dim,r_angle,axis)
      # Create VelocityLinear instance 
      Vel = genepy.VelocityLinear(Domain,u_norm,u_dir,u_type,u_angle,Rotation)

    Attributes
    ----------

    .. py:attribute:: norm
      :type: float
      :canonical: genepy.boundary_conditions.velocity.VelocityLinear.norm

      Velocity norm

    .. py:attribute:: alpha
      :type: float
      :canonical: genepy.boundary_conditions.velocity.VelocityLinear.alpha

      Velocity angle with the :math:`z` axis in radians
    
    .. py:attribute:: type
      :type: str
      :canonical: genepy.boundary_conditions.velocity.VelocityLinear.type

      Type of velocity field (``"extension"`` or ``"compression"``)

    .. py:attribute:: dir
      :type: int
      :canonical: genepy.boundary_conditions.velocity.VelocityLinear.dir

      Direction in which the velocity varies (``0``, ``1``, ``2``) defined from
      the given directions ``"x"``, ``"y"``, ``"z"`` respectively. 
    
    .. py:attribute:: uO
      :type: numpy.ndarray
      :canonical: genepy.boundary_conditions.velocity.VelocityLinear.uO

      Velocity vector :math:`\\mathbf u_O` at :math:`\\mathbf x = \\mathbf O`. Shape ``(dim,)``.
      Its value is computed at class initialization by the method :py:meth:`boundary_vector() <genepy.VelocityLinear.boundary_vector>`.

    .. py:attribute:: uL
      :type: numpy.ndarray
      :canonical: genepy.boundary_conditions.velocity.VelocityLinear.uL

      Velocity vector :math:`\\mathbf u_L` at :math:`\\mathbf x = \\mathbf L`. Shape ``(dim,)``.
      Its value is computed at class initialization by the method :py:meth:`boundary_vector() <genepy.VelocityLinear.boundary_vector>`.
    
    .. py:attribute:: a
      :type: numpy.ndarray
      :canonical: genepy.boundary_conditions.velocity.VelocityLinear.a

      Coefficients in vector form :math:`\\mathbf a` of the linear velocity function 
      :math:`\\mathbf u(x) = \\mathbf a x + \\mathbf b`. Shape ``(dim,)``.
      Its value is computed at class initialization by the method :py:meth:`velocity_coefficients() <genepy.VelocityLinear.velocity_coefficients>`.

    .. py:attribute:: b
      :type: numpy.ndarray
      :canonical: genepy.boundary_conditions.velocity.VelocityLinear.b

      Coefficients in vector form :math:`\\mathbf b` of the linear velocity function
      :math:`\\mathbf u(x) = \\mathbf a x + \\mathbf b`. Shape ``(dim,)``.
      Its value is computed at class initialization by the method :py:meth:`velocity_coefficients() <genepy.VelocityLinear.velocity_coefficients>`.

    .. py:attribute:: vertical_evaluated
      :type: bool

      Flag to check if the vertical velocity has been evaluated
    
    .. py:attribute:: u
      :type: numpy.ndarray
      :canonical: genepy.VelocityLinear.u

      Vector valued function of the velocity field. Shape ``(dim,)``. Contains sympy expressions.
    
    .. py:attribute:: grad_u
      :type: numpy.ndarray
      :canonical: genepy.VelocityLinear.grad_u

      Matrix of the gradient of the velocity field. Shape ``(dim,dim)``. Contains sympy expressions.
    
    .. py:attribute:: u_dir_horizontal
      :type: numpy.ndarray
      :canonical: genepy.VelocityLinear.u_dir_horizontal

      Horizontal vector of the velocity field orientation. Shape ``(dim,)``.
    
    .. py:attribute:: u_dir
      :type: numpy.ndarray
      :canonical: genepy.VelocityLinear.u_dir

      Full vector of the velocity field orientation. Shape ``(dim,)``.

    Methods
    -------
  """

  def __init__(self,Domain:domain.Domain,u_norm,variation_dir:str,velocity_type:str,u_angle:float=0.0,Rotation:rotation.Rotation=None) -> None:
    Velocity.__init__(self,Domain,Rotation)
    self.norm   = u_norm
    self.alpha  = u_angle
    self.type   = velocity_type

    if type(variation_dir) is not str:
      raise TypeError(f'variation_dir must be a string, found {type(variation_dir)}')
    
    dirmap = {"x":0,"y":1,"z":2}
    if variation_dir not in dirmap:
      raise ValueError(f"Invalid direction, possible values are \"x\", \"y\", \"z\", found: {variation_dir}")
    self.dir = dirmap[variation_dir]
    
    self.vertical_evaluated = False
    self.velocity_coefficients()
    self.u, self.grad_u   = self.evaluate_velocity_and_gradient_symbolic()
    self.u_dir_horizontal = self.get_velocity_orientation(horizontal=True,normalize=True)
    self.u_dir            = self.get_velocity_orientation(normalize=True)

  def __str__(self) -> str:
    s = f'{self.__class__.__name__}:\n'
    s += f'\tVelocity norm:        \t{self.norm}\n'
    s += f'\tVelocity angle:       \t{self.alpha} rad -> {np.rad2deg(self.alpha)} Deg\n'
    s += f'\tVelocity type:        \t{self.type}\n'
    s += f'\tVariation direction:  \t{self.sym_coor[self.dir]}\n'
    s += f'\tReferential rotation: \t{self.theta} rad -> {np.rad2deg(self.theta)}\n'
    s += f'\tRotation axis:        \t{self.axis}\n'
    return s

  def __repr__(self) -> str:
    """
    __repr__(self)
    Returns a string with the class name and attributes of the class.
    Can be used to print the class attributes.
    """
    attributes = vars(self)
    s  = f'{self.__class__.__name__}:\n'
    for attribute in attributes:
      if attribute == 'num_coor':
        s += f'\t{attribute}.shape:\t{attributes[attribute][0].shape}\n'
      else:
        s += f'\t{attribute}:\t{attributes[attribute]}\n'
    return s
  
  def boundary_vector(self,norm,alpha):
    """
    boundary_vector(self,norm,alpha)
    Computes a vector horizontal components from:

    - 1D & 2D: only the :attr:`norm <genepy.boundary_conditions.velocity.VelocityLinear.norm>` is used,
    - 3D: the :attr:`norm <genepy.boundary_conditions.velocity.VelocityLinear.norm>` 
      and angle :attr:`alpha <genepy.boundary_conditions.velocity.VelocityLinear.alpha>` 
      with the :math:`z` axis (North-South) such that:

    .. math::
      u_x &= \\sqrt{||\\mathbf u||^2 - u_z^2} \\\\
      u_y &= 0 \\\\
      u_z &= ||\\mathbf u|| \\cos(\\alpha)

    .. warning:: 
      These values are always positive, special attention is required 
      if different signs are needed. This is normally addressed with the 
      :attr:`type <genepy.boundary_conditions.velocity.VelocityLinear.type>` class attribute.

    :param float norm: velocity norm
    :param float alpha: angle with the :math:`z` axis in radians

    :return: boundary velocity, shape (dim,)
    :rtype: numpy.ndarray
    """
    if type(self.O[0]) is sp.Symbol: dtype = 'object'
    else:                            dtype = np.float64
    u = np.zeros(shape=(self.dim), dtype=dtype)
    if self.dim == 1 or self.dim == 2: 
      u[0] = norm
    elif self.dim == 3:
      u[2] = norm * np.cos(alpha)
      u[0] = np.sqrt(norm**2 - u[2]**2)
    return u
  
  def velocity_boundary(self):
    """
    velocity_boundary(self)
    Computes the boundary velocity vectors :py:attr:`uO <genepy.boundary_conditions.velocity.VelocityLinear.uO>` 
    and :py:attr:`uL <genepy.boundary_conditions.velocity.VelocityLinear.uL>` based on the 
    :py:attr:`norm <genepy.boundary_conditions.velocity.VelocityLinear.norm>` and the angle
    :py:attr:`alpha <genepy.boundary_conditions.velocity.VelocityLinear.alpha>` provided at class initialization.

    .. note::
      The velocity field is **symmetric** on the boundaries i.e., ``uL = -uO`` for **compression** and 
      ``uO = -uL`` for **extension**.

    """
    if self.type == "compression":
      self.uO = self.boundary_vector(self.norm,self.alpha)
      self.uL = -self.uO
    elif self.type == "extension":
      self.uL = self.boundary_vector(self.norm,self.alpha)
      self.uO = -self.uL
    else:
      raise RuntimeError(f'velocity_type can only be \"extension\" or \"compression\", found {self.type}')
    return

  def velocity_coefficients_1d(self,uO,uL,O,L):
    """
    velocity_coefficients_1d(self,uO,uL,O,L)
    computes velocity function coefficients :math:`a` and :math:`b` of the linear function 
    :math:`u(x) = ax + b` for a given direction :math:`x \\in [O,L]` such that:

    .. math::
      a &= \\frac{u_L - u_O}{L - O} \\\\
      b &= -aL + u_L
    
    with :math:`u_O` the velocity at :math:`x = O` 
    and :math:`u_L` the velocity at :math:`x = L` .

    :param float uO: velocity at the origin of the direction
    :param float uL: velocity at the max value of the direction
    :param float O: origin of the direction
    :param float L: max value of the direction

    :return: **a**: slope of the component of the velocity, 
             **b**: constant of the component of the velocity
    """
    a = (uL - uO) / (L - O)
    b = -a*L + uL
    return a,b
  
  def velocity_coefficients(self):
    """ 
    velocity_coefficients(self)
    For the vector valued velocity function :math:`\\mathbf u(x) = \\mathbf a x + \\mathbf b`,
    computes coefficients :math:`\\mathbf a` and :math:`\\mathbf b` such that:

    .. math:: 
      \\mathbf a &= \\frac{\\mathbf u_L - \\mathbf u_O}{\\mathbf L - \\mathbf O} \\\\
      \\mathbf b &= -\\mathbf a \\mathbf L + \\mathbf u_L
    
    with :math:`\\mathbf u_O` and :math:`\\mathbf u_L` computed by 
    :py:meth:`boundary_vector() <genepy.VelocityLinear.boundary_vector>` 
    and :math:`\\mathbf O` and :math:`\\mathbf L` the origin and max values of the domain.
    Each component of the velocity field is computed by calling 
    :py:meth:`velocity_coefficients_1d() <genepy.VelocityLinear.velocity_coefficients_1d>`.
    """
    if type(self.O[0]) is sp.Symbol: dtype = 'object'
    else:                            dtype = np.float64

    self.a = np.zeros(shape=(self.dim), dtype=dtype)
    self.b = np.zeros(shape=(self.dim), dtype=dtype)
    
    self.velocity_boundary()

    for d in range(self.dim):
      self.a[d], self.b[d] = self.velocity_coefficients_1d(self.uO[d],self.uL[d],self.O[self.dir],self.L[self.dir])
    return
  
  def velocity_function_1d(self,x,a,b):
    """
    velocity_function_1d(self,x,a,b)
    Computes the linear velocity field in 1 direction (scalar valued function) 
    such that :math:`u(x) = ax + b` for the given direction :math:`x` and 
    coefficients :math:`a` and :math:`b` computed by the method 
    :py:meth:`velocity_coefficients_1d <genepy.VelocityLinear.velocity_coefficients_1d>`.

    :param x: coordinate of the direction in which the velocity varies
    :param a: slope of the component of the velocity
    :param b: constant of the component of the velocity

    :return: scalar valued function of the velocity field
    """
    return a*x + b
  
  def velocity_function(self,x):
    """
    velocity_function(self,x)
    computes the vector valued linear velocity function
    such that :math:`\\mathbf u(x) = \\mathbf a x + \\mathbf b`.
    Calls the method :meth:`velocity_function_1d` for each direction.

    :param x: coordinates of the direction in which the velocity varies

    :return: vector valued velocity function
    """
    u = []
    for d in range(self.dim):
      u.append(self.velocity_function_1d(x,self.a[d],self.b[d]))
    return u
  
  def evaluate_velocity_symbolic(self):
    """
    evaluate_velocity_symbolic(self)
    Evaluates the velocity field in symbolic form returning 
    a vector valued function of the form :math:`\\mathbf u (x,y,z) = \\mathbf a x + \\mathbf b y + \\mathbf c z + \\mathbf d`
    and shape ``(1,dim)``.

    If a rotation is required i.e., an instance of the :py:class:`Rotation <genepy.Rotation>`
    class with a non-zero angle is provided, the velocity is evaluated as:

    .. math::
      \\mathbf u_R (\\mathbf x) = \\boldsymbol R \\mathbf u (\\mathbf x_R)

    where :math:`\\mathbf x` is the non-rotated coordinate system, 
    :math:`\\boldsymbol R` is the rotation matrix,
    :math:`\\mathbf x_R` is the rotated coordinate system,
    :math:`\\mathbf u` is the velocity field before rotation and
    :math:`\\mathbf u_R` is the rotated velocity field. 

    The rotation of the coordinate system is performed by the method 
    :meth:`rotate_referential() <genepy.Rotation.rotate_referential>`.
    The evaluation of the velocity field with the (rotated) coordinate system is done with
    the method :meth:`velocity_function() <genepy.VelocityLinear.velocity_function>`.
    """
    # create a numpy array of shape (1,dim) with the symbolic coordinates
    coor = np.array([[*self.sym_coor]], dtype='object')
    # rotate referential
    coor_R = self.rotate_referential(coor,self.O,self.L,ccw=False)
    # evaluate velocity with the rotated coordinate
    u = np.zeros(shape=(1,self.dim), dtype='object')
    u[0,:] = self.velocity_function(coor_R[0,self.dir])
    # rotate the velocity field
    R   = self.rotation_matrix()
    u_R = self.rotate_vector(R,u,ccw=True)
    # reshape for a 1D array
    u_R.shape = (self.dim,)
    return u_R
  
  def evaluate_velocity_and_gradient_symbolic(self):
    """
    evaluate_velocity_and_gradient_symbolic(self)
    Calls the methods:
    
    - :py:meth:`evaluate_velocity_symbolic <genepy.VelocityLinear.evaluate_velocity_symbolic>` 
      to evaluate the horizontal components of the velocity 
    - :py:meth:`evaluate_vertical_velocity <genepy.VelocityLinear.evaluate_vertical_velocity>` 
      to evaluate the vertical component of the velocity
    - :py:meth:`evaluate_gradient <genepy.Velocity.evaluate_gradient>` 
      to evaluate the gradient of the velocity vector  
    
    in symbolic form.
    """
    u      = self.evaluate_velocity_symbolic()
    u[1]   = self.evaluate_vertical_velocity(self.sym_coor[1],u=u)
    grad_u = self.evaluate_gradient(u)
    return u,grad_u
  
  def evaluate_vertical_velocity_coefficients(self,u=None):
    """
    Evaluates the coefficients of the vertical component of the velocity function
    :math:`a` and :math:`b` such that :math:`u_y(y) = a y + b`. 
    To evaluate the coefficients, the following boundary conditions are used

    .. math::

      u_y(L_y) &= 0 \\\\
      u_y(O_y) &= \\frac{1}{S_{xz}} \\sum_f \\int_S \\mathbf u \\cdot \\mathbf n \\, dS \\\\

    where :math:`S_{xz}` is the surface of the bottom face of the domain, 
    :math:`O_y` and :math:`L_y` are the minimum and maximum values 
    of the domain in the :math:`y` direction. 

    :math:`a` and :math:`b` are then computed with 
    :meth:`velocity_coefficients_1d() <genepy.VelocityLinear.velocity_coefficients_1d>` 

    :param u: **(Optional)**, vector valued function of the velocity field.
              If not provided, the symbolic velocity field is evaluated with :meth:`evaluate_velocity_symbolic`
    
    :return: **a**: slope of the vertical component of the velocity,
             **b**: constant of the vertical component of the velocity
    """
    if u is None:
      u = self.evaluate_velocity_symbolic()
    int_u_dot_n = self.evaluate_int_u_dot_n_faces(u)
    iudn = 0.0
    for face in int_u_dot_n:
      iudn += int_u_dot_n[face]
    bottom_surface = (self.L[0] - self.O[0])
    if self.dim == 3: bottom_surface *= (self.L[2] - self.O[2])
    self.uO[1] = iudn / bottom_surface
    self.uL[1] = 0.0

    self.a[1],self.b[1] = self.velocity_coefficients_1d(self.uO[1],self.uL[1],self.O[1],self.L[1])
    self.vertical_evaluated = True
    return

  def evaluate_vertical_velocity(self,y,u=None):
    """
    Evaluates the vertical component of the velocity field.
    Can be used for numerical or symbolic evaluation.
    Calls the methods:

    - :py:meth:`evaluate_vertical_velocity_coefficients() <genepy.VelocityLinear.evaluate_vertical_velocity_coefficients>`
    - :py:meth:`velocity_function_1d() <genepy.VelocityLinear.velocity_function_1d>`

    :param y: coordinate of the vertical direction
    :param u: **(Optional)**, vector valued function of the velocity field.
              If not provided, the symbolic velocity field is evaluated with :meth:`evaluate_velocity_symbolic`

    :return: vertical component of the velocity field
    """
    if self.vertical_evaluated == False:
      self.evaluate_vertical_velocity_coefficients(u=u)
    uy = self.velocity_function_1d(y,self.a[1],self.b[1])
    return uy
  
  def evaluate_velocity_numeric(self):
    """
    evaluate_velocity_numeric(self)
    Evaluates the velocity field numerically
    """
    coor = self.shape_coor()
    # rotate referential
    coor_R = self.rotate_referential(coor,self.O,self.L,ccw=False)
    # evaluate velocity with the rotated coordinate
    u = self.velocity_function(coor_R[:,self.dir])
    # convert to numpy array of shape (npoints, dim)
    u = np.asarray(u, dtype=np.float64).T
    # rotate the velocity field
    R = self.rotation_matrix()
    u_R = self.rotate_vector(R,u,ccw=True)
    u_R[:,1] = self.evaluate_vertical_velocity(coor[:,1])
    return u_R
  
  def plot_velocity_matplotlib(self):
    u = self.evaluate_velocity_numeric()
    X,Z = self.num_coor[0],self.num_coor[self.dim-1]
    ux = u[:,0].reshape(self.n)
    uz = u[:,self.dim-1].reshape(self.n)

    u_norm = np.sqrt(ux**2 + uz**2)

    _,ax = plt.subplots()
    ax.contourf(X,Z,u_norm,100,cmap='Blues')
    ax.quiver(X,Z,ux,uz)
    ax.axis('equal')
    plt.show()
    return
  
class VelocityLinearAsymmetric(VelocityLinear):
  """
  .. py:class:: VelocityLinearAsymmetric(Domain, u_normL, u_normO, variation_dir, velocity_type, u_angleL=0.0, u_angleO=0.0, Rotation=None)

    Inherits from :py:class:`VelocityLinear <genepy.VelocityLinear>` 
    to evaluate a linear velocity field with different velocities at the boundaries.
    Both the norm and the angle can be different at the boundaries.

    :param Domain Domain: domain in which the velocity field is evaluated
    :param float u_normL: velocity norm of the vector at the boundary of minimum coordinate
    :param float u_normO: velocity norm of the vector at the boundary of maximum coordinate
    :param str variation_dir: direction in which the velocity varies (``"x"``, ``"y"``, ``"z"``)
    :param str velocity_type: velocity field orientation, 
                              ``"extension"`` velocity is directed outward the domain, 
                              ``"compression"`` velocity is directed inward the domain
    :param float u_angleL: **(Optional)**  angle in radians of the velocity field at the boundary of minimum coordinate
    :param float u_angleO: **(Optional)**  angle in radians of the velocity field at the boundary of maximum coordinate
    :param Rotation Rotation: **(Optional)** Rotation class instance to rotate the referential

    Example:
    --------

    For an asymetric orthogonal shortening in the :math:`z` direction from 1 cm/a to 0 cm/a at the boundaries,
    assuming that the class :class:`Domain <genepy.Domain>` 
    has been instanciated as ``domain``:

    .. code-block:: python

      cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0) # cm/a to m/s conversion
      u_normL = 1.0*cma2ms
      u_normO = 0.0
      u_dir   = "z"
      u_type  = "compression"
      Vel = gp.VelocityLinearAsymmetric(domain,u_normL,u_normO,u_dir,u_type)

      
    Attributes
    ----------

    .. py:attribute:: normL
      :type: float
      :canonical: genepy.boundary_conditions.velocity.VelocityLinearAsymmetric.normL

      Velocity norm at the boundary of minimum coordinate

    .. py:attribute:: normO
      :type: float
      :canonical: genepy.boundary_conditions.velocity.VelocityLinearAsymmetric.normO

      Velocity norm at the boundary of maximum coordinate

    .. py:attribute:: alphaL
      :type: float
      :canonical: genepy.boundary_conditions.velocity.VelocityLinearAsymmetric.alphaL

      Velocity angle with the :math:`z` axis in radians at the boundary of minimum coordinate. Default: ``0.0``.

    .. py:attribute:: alphaO
      :type: float
      :canonical: genepy.boundary_conditions.velocity.VelocityLinearAsymmetric.alphaO

      Velocity angle with the :math:`z` axis in radians at the boundary of maximum coordinate. Default: ``0.0``.

  """
  def __init__(self, Domain: domain.Domain, 
               u_normL:float,u_normO:float,variation_dir:str,velocity_type:str,u_angleL:float=0.0,u_angleO:float=0.0,
               Rotation: rotation.Rotation = None) -> None:
    Velocity.__init__(self,Domain,Rotation)
    self.normL  = u_normL
    self.normO  = u_normO
    self.alphaL = u_angleL
    self.alphaO = u_angleO
    self.type   = velocity_type

    if type(variation_dir) is not str:
      raise TypeError(f'variation_dir must be a string, found {type(variation_dir)}')
    
    dirmap = {"x":0,"y":1,"z":2}
    if variation_dir not in dirmap:
      raise ValueError(f"Invalid direction, possible values are \"x\", \"y\", \"z\", found: {variation_dir}")
    self.dir = dirmap[variation_dir]

    self.vertical_evaluated = False
    self.velocity_coefficients()
    self.u, self.grad_u   = self.evaluate_velocity_and_gradient_symbolic()
    self.u_dir_horizontal = self.get_velocity_orientation(horizontal=True,normalize=True)
    self.u_dir            = self.get_velocity_orientation(normalize=True)
    return

  def velocity_boundary(self):
    """
    velocity_boundary(self)
    Computes the boundary velocity vectors :py:attr:`uO <genepy.boundary_conditions.velocity.VelocityLinear.uO>` 
    and :py:attr:`uL <genepy.boundary_conditions.velocity.VelocityLinear.uL>` based on the 
    :py:attr:`normL <genepy.boundary_conditions.velocity.VelocityLinearAsymmetric.normL>`,
    :py:attr:`normO <genepy.boundary_conditions.velocity.VelocityLinearAsymmetric.normO>` and the angles
    :py:attr:`alphaL <genepy.boundary_conditions.velocity.VelocityLinearAsymmetric.alphaL>`,
    :py:attr:`alphaO <genepy.boundary_conditions.velocity.VelocityLinearAsymmetric.alphaO>` provided at class initialization.

    .. note::
      In compression :py:attr:`uL <genepy.boundary_conditions.velocity.VelocityLinear.uL>` 
      is directed inward the domain i.e., negative, and in extension 
      :py:attr:`uO <genepy.boundary_conditions.velocity.VelocityLinear.uO>` 
      is directed outward the domain i.e., positive.
      
    """
    if self.type == "compression":
      self.uO =  self.boundary_vector(self.normO,self.alphaO)
      self.uL = -self.boundary_vector(self.normL,self.alphaL)
    elif self.type == "extension":
      self.uL =  self.boundary_vector(self.normL,self.alphaL)
      self.uO = -self.boundary_vector(self.normO,self.alphaO)
    else:
      raise RuntimeError(f'velocity_type can only be \"extension\" or \"compression\", found {self.type}')
    return
  
class VelocityTimeDependant(Velocity):
  """
  .. py:class:: VelocityTimeDependant(Domain, Rotation=None)

    Class to construct and evaluate time dependant velocity function.
    This class is the parent class of all time dependant velocity classes 
    and inherits from the class :py:class:`Velocity <genepy.Velocity>`.

    :param Domain Domain: domain in which the velocity field is evaluated
    :param Rotation Rotation: **(Optional)** Rotation class instance to rotate the referential

    Methods
    -------
  """
  def __init__(self, Domain:domain.Domain, Rotation:rotation.Rotation=None) -> None:
    Velocity.__init__(self, Domain, Rotation)
    self.time_sym = sp.symbols('t')
    return

  def sum_functions(self,user_func:list,user_args:list):
    """
    sum_functions(self,user_func,user_args)
    Abstract method to sum user defined functions. 
    The functions to sum should return objects of the same shape and type.

    Considering the functions :math:`f_1(x), f_2(x), \\ldots , f_n(x)` this method computes:
    
    .. math:: 
      u(x) = \\sum_{i=1}^{n} f_i(x)
    
    where :math:`x` represents the variable(s) of the function.

    :param list user_func: list of user defined functions
    :param list user_args: list of arguments to pass to the functions
    """
    u = user_func[0](*(user_args[0]))
    nf = len(user_func)
    if nf > 1:
      for n in range(1,nf):
        u += user_func[n](*(user_args[n]))
    return u

class VelocityInversion(VelocityTimeDependant):
  """
  .. py:class:: VelocityInversion(Domain, phase1, phase2, breakpoints, slopes, Rotation=None)

    Class to construct and evaluate a time dependant velocity function with two phases by summing
    two :py:meth:`arctangent functions <genepy.VelocityInversion.velocity_function>`.
    
    :param Domain Domain: domain in which the velocity field is evaluated
    :param VelocityLinear phase1: :py:class:`Velocity <genepy.VelocityLinear>` class instance of the first phase
    :param VelocityLinear phase2: :py:class:`Velocity <genepy.VelocityLinear>` class instance of the second phase
    :param list breakpoints: list of two breakpoints in time
    :param list slopes: list of two slopes at the breakpoints
    :param Rotation Rotation: **(Optional)** Rotation class instance to rotate the referential

    
    Attributes
    ----------

    .. py:attribute:: phases
      :type: list

      List containing the 2 class instances of the velocity phases. Shape ``[phase1:Velocity, phase2:Velocity]``.
    
    .. py:attribute:: breakpoints
      :type: list

      List of the two breakpoints in time. Shape ``[breakpoint1:float, breakpoint2:float]``
    
    .. py:attribute:: slopes
      :type: list

      List of the two slopes at the breakpoints. Shape ``[slope1:float, slope2:float]``

    .. py:attribute:: u
      :type: numpy.ndarray

      Vector valued function of the velocity field. Shape ``(dim,)``. Contains sympy expressions.
    
    .. py:attribute:: grad_u
      :type: numpy.ndarray

      Matrix of the gradient of the velocity field. Shape ``(dim,dim)``. Contains sympy expressions.
    
    .. py:attribute:: u_dir_horizontal
      :type: list

      List of horizontal vectors of the velocity field orientation for each phase. Shape ``[(dim,), (dim,)]``.
    
    .. py:attribute:: u_dir
      :type: numpy.ndarray

      List of the full vector of the velocity field orientation for each phase. Shape ``[(dim,), (dim,)]``.

    Methods
    -------
  """
  def __init__(self, Domain:domain.Domain, 
               phase1:VelocityLinear, 
               phase2:VelocityLinear, 
               breakpoints, 
               slopes,
               Rotation:rotation.Rotation=None) -> None:
    VelocityTimeDependant.__init__(self, Domain, Rotation)
    self.phases      = [phase1, phase2]
    self.breakpoints = breakpoints
    self.slopes      = slopes

    self.u, self.grad_u   = self.evaluate_velocity_and_gradient_symbolic()
    self.u_dir_horizontal = [self.phases[0].u_dir_horizontal, self.phases[1].u_dir_horizontal]
    self.u_dir            = [self.phases[0].u_dir, self.phases[1].u_dir]
    return

  def velocity_function(self, time, bound, **kwargs):
    """
    velocity_function(self, time, bound)
    Evaluates a time dependant velocity function :math:`u(t)` with two phases by summing two arctangent functions such that

    .. math::
      f_1(t) &= b_1 \\left( \\frac{1}{2} - \\frac{\\tan^{-1} \\left(s_1(t-t_1) \\right)}{\\pi} \\right) \\\\
      f_2(t) &= b_2 \\left( \\frac{1}{2} + \\frac{\\tan^{-1} \\left(s_2(t-t_2) \\right)}{\\pi} \\right) \\\\
      u(t) &= f_1(t) + f_2(t)

    where :math:`t_i` are the breakpoints in time, :math:`s_i` are the arctangent slopes at the breakpoints 
    and :math:`b_i` are the bounds of the functions (min and max values i.e., the steady state velocities of each phase).

    :param float time: time at which the velocity function is evaluated, can be symbolic or numeric
    :param list bound: list of two bounds of the velocity function

    :return: time-dependant velocity function
    """
    def phase_1(t, t0, s, b):
      fac = 2.0/np.pi
      func = 0.5*( -fac*b*sp.atan(s*(t-t0)) + b )
      return func
    
    def phase_2(t, t0, s, b):
      fac = 2.0/np.pi
      func = 0.5*( fac*b*sp.atan(s*(t-t0)) + b )
      return func
    
    user_func = [phase_1, phase_2]
    user_args = [(time, self.breakpoints[0], self.slopes[0], bound[0]), (time, self.breakpoints[1], self.slopes[1], bound[1])]
    u_t = self.sum_functions(user_func, user_args)
    return u_t
  
  def velocity_function_derivative(self,time,bound,**kwargs):
    """
    velocity_function_derivative(self,time,bound)
    Evaluates the derivative with respect to time of the velocity function.
    The derivative is computed by summing the derivatives of the two arctangent functions defined in 
    :py:meth:`velocity_function <genepy.VelocityInversion.velocity_function>` such that

    .. math:: 
      f'_1(t) &= -\\frac{b_1 s_1}{\\pi (1 + (s_1(t-t_1))^2)} \\\\
      f'_2(t) &=  \\frac{b_2 s_2}{\\pi (1 + (s_2(t-t_2))^2)}

    :param float time: time at which the function is evaluated, can be symbolic or numeric
    :param list bound: list of two bounds of the stady-state velocity functions

    :return: derivative with respect to time of the velocity function
    """
    def phase_1_derivative(t, t0, s, b):
      fac = 1.0/np.pi
      func = -fac*b*s/(1.0 + (s*(t-t0))**2)
      return func

    def phase_2_derivative(t, t0, s, b):
      fac = 1.0/np.pi
      func = fac*b*s/(1.0 + (s*(t-t0))**2)
      return func
    
    user_func = [phase_1_derivative, phase_2_derivative]
    user_args = [(time, self.breakpoints[0], self.slopes[0], bound[0]), 
                 (time, self.breakpoints[1], self.slopes[1], bound[1])]
    du_dt = self.sum_functions(user_func, user_args)
    return du_dt

  def evaluate_velocity_symbolic(self):
    """
    evaluate_velocity_symbolic(self)
    Evaluates the time dependant velocity function in symbolic form.
    Calls the methods :py:meth:`evaluate_velocity_symbolic <genepy.VelocityLinear.evaluate_velocity_symbolic>` 
    and :py:meth:`evaluate_vertical_velocity <genepy.VelocityLinear.evaluate_vertical_velocity>` for each phase.
    Then evaluates the time dependant velocity function with the method 
    :py:meth:`velocity_function <genepy.VelocityInversion.velocity_function>`.

    :return: time dependant velocity function
    """
    u1 = self.phases[0].u
    u2 = self.phases[1].u

    u_t = []
    for d in range(self.phases[0].dim):
      u_bounds = [ u1[d], u2[d] ]
      u_t.append(self.velocity_function(self.time_sym,u_bounds))
    return u_t
  
  def evaluate_velocity_and_gradient_symbolic(self):
    """
    evaluate_velocity_and_gradient_symbolic(self)
    Calls the methods:

    - :meth:`evaluate_velocity_symbolic <genepy.VelocityInversion.evaluate_velocity_symbolic>` 
      to evaluate the time dependant velocity function
    - :meth:`evaluate_gradient <genepy.Velocity.evaluate_gradient>` 
      to evaluate the gradient of the velocity function

    :return: ``u, grad_u``, time dependant velocity function and its gradient
    """
    u = self.evaluate_velocity_symbolic()
    grad_u = self.evaluate_gradient(u)
    return u,grad_u

  def get_time_zero_velocity(self,report=False):
    """
    get_time_zero_velocity(self,report=False)
    Computes the time at which the velocity function evaluates to zero.
    The result is obtained by solving the equation :math:`f(t) = 0` with the iterative Newton-Raphson method.
    The initial guess is set to the first breakpoint in time because it is the most likely to converge (steep part of the function).

    :param bool report: if True, the Newton-Raphson method reports the convergence of the solution (default: ``False``)

    :return: time at which the velocity function evaluates to zero
    """
    # evaluate 1D function
    u_1 = self.phases[0].norm
    if self.phases[0].type == "compression": u_1 = -u_1
    u_2 = self.phases[1].norm
    if self.phases[1].type == "compression": u_2 = -u_2

    u_bounds = [u_1, u_2]
    # scale for numerical stability
    u_scale = 10**(-np.floor(np.log10(u_1)))
    # use one of the breakpoint as initial guess, it should converge fast
    t_guess  = self.breakpoints[0]
    t0 = newton_raphson(self.velocity_function,
                        self.velocity_function_derivative,
                        t_guess,tol=1.0e-6,max_iter=10,scaling=u_scale,report=report,
                        bound=u_bounds)
    return t0

  def plot_1D_velocity(self,time):
    """
    plot_1D_velocity(self,time)
    Plots the velocity function in 1D over time.
    If the spatial dimension is more than 1 (i.e., 2D and 3D)
    only the norm of the velocity function is evaluated to allow plotting.
    The sign convention for the plot is:

    - extension: positive velocity
    - compression: negative velocity

    Displays two plots:

    - the velocity over time in m/s
    - the velocity over time in cm/a

    On both plots the point at which the velocity evaluates to zero is marked with a red circle.

    :param numpy.ndarray time: time array
    """
    # evaluate 1D function
    u_1 = self.phases[0].norm
    if self.phases[0].type == "compression": u_1 = -u_1
    u_2 = self.phases[1].norm
    if self.phases[1].type == "compression": u_2 = -u_2

    u_bounds = [u_1, u_2]
    u_t = np.zeros(shape=(len(time)), dtype=np.float64)
    for i,t in enumerate(time):
      u_t[i] = self.velocity_function(t,u_bounds)
    
    t0 = self.get_time_zero_velocity()

    _,ax = plt.subplots(1,2)
    ax[0].plot(time,u_t,'tab:blue')
    ax[0].plot(t0,0.0,'ro',markeredgecolor='black')
    ax[0].set_xlabel('Time (s)')
    ax[0].set_ylabel('Velocity (m/s)')
    ax[0].set_title('Velocity in m/s over time in s')

    Myr2sec = (3600.0 * 24.0 * 365.0) * 1e6
    cma2mps = 1e-2 / (3600.0 * 24.0 * 365.0)
    ax[1].plot(time/Myr2sec,u_t/cma2mps,'tab:green')
    ax[1].plot(t0/Myr2sec,0.0,'ro',markeredgecolor='black')
    ax[1].set_xlabel('Time (Myr)')
    ax[1].set_ylabel('Velocity (cm/a)')
    ax[1].set_title('Velocity in cm/a over time in Myr')

    plt.show()
    return
  
  def paraview_velocity_inversion(self,writer:writers.WriteVTS,time,root:str,pvd:str):
    import os
    u1 = self.phases[0].evaluate_velocity_numeric()
    u2 = self.phases[1].evaluate_velocity_numeric()

    u_bounds = []
    for d in range(self.phases[0].dim):
      u_bounds.append([ u1[:,d], u2[:,d] ])

    step = 0
    writer.pvd_fname = os.path.join(root,pvd)
    writer.open_pvd()
    for t in time:
      output = f"velocity_{step}.vts"
      writer.vtk_fname = os.path.join(root,output)

      u_t = []
      for d in range(self.phases[0].dim):
        u_t.append(self.velocity_function(t,u_bounds[d]))
      u_t = np.asarray(u_t, dtype=np.float64).T
      writer.point_data = { "u": u_t }
      writer.write_vts()
      writer.append_pvd(t,output)
      step += 1
    writer.close_pvd()                                          
    return

class VelocityCompose(Velocity):
  """
  .. py:class:: VelocityCompose(Domain, Velocities, Rotation=None, method="sum")

    Class to construct and evaluate a velocity function composed by multiple velocity functions.
    The composition method is defined by the ``method`` parameter.

    :param Domain Domain: domain in which the velocity field is evaluated
    :param list Velocities: list of :py:class:`Velocity <genepy.Velocity>` class instances
    :param Rotation Rotation: **(Optional)** Rotation
    :param str method: **(Optional)** method to compose the velocity functions, default is ``"sum"``

    Attributes
    ----------

    .. py:attribute:: nfuncs
      :type: int
      :canonical: genepy.boundary_conditions.velocity.VelocityCompose.nfuncs

      Number of velocity functions to compose

    .. py:attribute:: velocities
      :type: list
      :canonical: genepy.boundary_conditions.velocity.VelocityCompose.velocities

      List of :py:class:`velocity <genepy.Velocity>` class instances of velocity functions to compose.

    .. py:attribute:: method
      :type: str
      :canonical: genepy.boundary_conditions.velocity.VelocityCompose.method

      Method to compose the velocity functions. Default is ``"sum"``.

    Methods
    -------
  """
  def __init__(self, Domain:domain.Domain, Velocities:list[Velocity], Rotation:rotation.Rotation=None, method="sum") -> None:
    Velocity.__init__(self, Domain, Rotation)
    self.nfuncs     = len(Velocities)
    self.velocities = Velocities
    if type(method) is not str: raise TypeError(f'method must be a string, found {type(method)}')
    self.method = method
    self.u = np.zeros(shape=(self.dim), dtype='object')

    if self.method == "sum": self.sum_velocities()
    else: raise ValueError(f'Invalid method, possible values are \"sum\", found {method}')
    
    self.grad_u           = self.evaluate_gradient(self.u)
    self.u_dir_horizontal = self.get_velocity_orientation(horizontal=True,normalize=True)
    self.u_dir            = self.get_velocity_orientation(normalize=True)
    return

  def sum_velocities(self):
    """
    sum_velocities(self)
    Sums the velocity functions of :py:class:`velocity <genepy.Velocity>` class instances from the list 
    :py:attr:`velocities <genepy.boundary_conditions.velocity.VelocityCompose.velocities>`.
    """
    self.u = sum(v.u for v in self.velocities)
    return