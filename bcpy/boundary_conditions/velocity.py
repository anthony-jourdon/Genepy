from bcpy.initial_conditions import domain
from bcpy import rotation
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

class Velocity(domain.Domain,rotation.Rotation):
  """
  .. py:class:: Velocity(Domain, u_norm, variation_dir, velocity_type, u_angle=0.0, Rotation=None)

    Class to evaluate symbolic and numeric linear velocity function and its derivatives in space.
    The velocity field is defined as a linear function of the coordinates. 
    It can be used for 2D and 3D domains.
    Given a horizontal velocity field, the vertical velocity is computed by integrating the dot product of the velocity field with the normal vector over the faces of the domain.
    The class inherits from :py:class:`Domain <bcpy.initial_conditions.domain.Domain>` and :py:class:`Rotation <bcpy.rotation.Rotation>`

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
    Assuming the class :class:`Domain <bcpy.initial_conditions.domain.Domain>` as been instanciated as ``Domain``
    
    **Without rotation** of the referential

    .. code-block:: python

      u_norm  = 1.0              # horizontal velocity norm
      u_angle = np.deg2rad(45.0) # velocity angle \in [-pi/2, pi/2]. If 0, can be ommited
      u_dir   = "z"              # direction in which velocity varies
      u_type  = "extension"      # velocity orientation
      # Create Velocity instance 
      Vel = bcpy.Velocity(Domain,u_norm,u_dir,u_type,u_angle)
    
    **With rotation** of the referential

    .. code-block::  python

      # Rotation of the referential
      r_angle = np.deg2rad(-15.0)                   # Rotation angle
      axis    = np.array([0,1,0], dtype=np.float64) # Rotation axis
      # Create Rotation instance
      Rotation = bp.Rotation(dim,r_angle,axis)
      # Create Velocity instance 
      Vel = bcpy.Velocity(Domain,u_norm,u_dir,u_type,u_angle)

    Attributes
    ----------

    .. py:attribute:: norm
      :type: float
      :canonical: bcpy.boundary_conditions.velocity.Velocity.norm

      Velocity norm

    .. py:attribute:: alpha
      :type: float

      Velocity angle with the z axis in radians
    
    .. py:attribute:: type
      :type: str

      Type of velocity field (``"extension"``, ``"compression"``)

    .. py:attribute:: dir
      :type: int

      Direction in which the velocity varies (0, 1, 2)
    
    .. py:attribute:: uO
      :type: numpy.ndarray
      :canonical: bcpy.boundary_conditions.velocity.Velocity.uO

      Velocity vector :math:`\\mathbf u_O` at :math:`\\mathbf x = \\mathbf O`. Shape ``(dim,)``.
      Its value is computed at class initialization by the method :py:meth:`boundary_vector() <bcpy.Velocity.boundary_vector>`.

    .. py:attribute:: uL
      :type: numpy.ndarray
      :canonical: bcpy.boundary_conditions.velocity.Velocity.uL

      Velocity vector :math:`\\mathbf u_L` at :math:`\\mathbf x = \\mathbf L`. Shape ``(dim,)``.
      Its value is computed at class initialization by the method :py:meth:`boundary_vector() <bcpy.Velocity.boundary_vector>`.
    
    .. py:attribute:: a
      :type: numpy.ndarray
      :canonical: bcpy.boundary_conditions.velocity.Velocity.a

      Coefficients in vector form :math:`\\mathbf a` of the linear velocity function 
      :math:`\\mathbf u(x) = \\mathbf a x + \\mathbf b`. Shape ``(dim,)``.
      Its value is computed at class initialization by the method :py:meth:`velocity_coefficients() <bcpy.Velocity.velocity_coefficients>`.

    .. py:attribute:: b
      :type: numpy.ndarray
      :canonical: bcpy.boundary_conditions.velocity.Velocity.b

      Coefficients in vector form :math:`\\mathbf b` of the linear velocity function
      :math:`\\mathbf u(x) = \\mathbf a x + \\mathbf b`. Shape ``(dim,)``.
      Its value is computed at class initialization by the method :py:meth:`velocity_coefficients() <bcpy.Velocity.velocity_coefficients>`.

    .. py:attribute:: vertical_evaluated
      :type: bool

      Flag to check if the vertical velocity has been evaluated
    
    Methods
    -------

  """

  def __init__(self,Domain:domain.Domain,u_norm,variation_dir:str,velocity_type:str,u_angle:float=0.0,Rotation:rotation.Rotation=None) -> None:
    self.norm   = u_norm
    self.alpha  = u_angle
    self.type   = velocity_type

    if type(variation_dir) is not str:
      raise TypeError(f'variation_dir must be a string, found {type(variation_dir)}')
    
    for d in variation_dir:
      if   d == "x": self.dir = 0
      elif d == "y": self.dir = 1
      elif d == "z": self.dir = 2
      else: raise ValueError(f"Invalid direction, possible values are \"x\", \"y\", \"z\", found: {d}")

    domain.Domain.__init__(self,Domain.dim,Domain.O,Domain.L,Domain.n,coor=Domain.num_coor)
    # set rotation angle to zero if Rotation class is not provided
    if Rotation is None: rotation.Rotation.__init__(self,Domain.dim,0.0,np.array([0,1,0]))
    else:                rotation.Rotation.__init__(self,Domain.dim,Rotation.theta,Rotation.axis)

    self.vertical_evaluated = False
    self.velocity_coefficients()

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
  
  def report_symbolic_functions(self,u,grad_u,uL) -> str:
    """
    report_symbolic_functions(self,u,grad_u,uL)
    Returns a string with the symbolic velocity function, 
    its derivatives and the boundary velocity orientation.
    Can be used with print() or written to a file.

    :param numpy.ndarray u: sympy vector valued function of the velocity field
    :param numpy.ndarray grad_u: sympy matrix of the derivatives of the velocity field shape (dim,dim)
    :param numpy.ndarray uL: sympy orientation vector of the velocity field at the boundary

    :return: string with the symbolic velocity function, its derivatives and the boundary velocity orientation.
    :rtype: str
    """
    s = f"Symbolic velocity function:\n"
    for i in range(self.dim):
      s += f"\tu{self.sym_coor[i]}{self.sym_coor} = {u[0,i]}\n"
    s += f"Derivatives of the velocity function:\n"
    for i in range(self.dim):
      for j in range(self.dim):
        s += f"\tdu{self.sym_coor[i]}/d{self.sym_coor[j]} = "+str(grad_u[i,j])+"\n"
    s += f"Boundary velocity orientation:\n"
    for i in range(self.dim):
      s += f"\tuL{self.sym_coor[i]} = {uL[i]}\n"
    return s
  
  def boundary_vector(self):
    """
    boundary_vector(self)
    Computes the vector horizontal components.
    1D & 2D: only the norm is used,
    3D:      based on the norm and angle with z axis (North-South) such that:

    .. math::
      u_x &= \\sqrt{||\\mathbf u||^2 - u_z^2} \\\\
      u_y &= 0 \\\\
      u_z &= ||\\mathbf u|| \\cos(\\alpha)

    .. warning:: 
      These values are always positive, special attention is required 
      if different signs are needed.

    :return: boundary velocity, shape (dim,)
    :rtype: numpy.ndarray
    """
    u = np.zeros(shape=(self.dim), dtype=np.float64)
    if self.dim == 1 or self.dim == 2: 
      u[0] = self.norm
    elif self.dim == 3:
      u[2] = self.norm * np.cos(self.alpha)
      u[0] = np.sqrt(self.norm**2 - u[2]**2)
    return u
  
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
    :py:meth:`boundary_vector() <bcpy.boundary_conditions.velocity.Velocity.boundary_vector>` 
    and :math:`\\mathbf O` and :math:`\\mathbf L` the origin and max values of the domain.
    Each component of the velocity field is computed by calling 
    :py:meth:`velocity_coefficients_1d() <bcpy.boundary_conditions.velocity.Velocity.velocity_coefficients_1d>`.
    """
    self.a = np.zeros(shape=(self.dim), dtype=np.float64)
    self.b = np.zeros(shape=(self.dim), dtype=np.float64)

    if self.type == "compression":
      self.uO = self.boundary_vector()
      self.uL = -self.uO
    elif self.type == "extension":
      self.uL = self.boundary_vector()
      self.uO = -self.uL
    else:
      raise RuntimeError(f'velocity_type can only be \"extension\" or \"compression\", found {self.type}')
    
    for d in range(self.dim):
      self.a[d], self.b[d] = self.velocity_coefficients_1d(self.uO[d],self.uL[d],self.O[self.dir],self.L[self.dir])
    return
  
  def linear_velocity_1d(self,x,a,b):
    """
    linear_velocity_1d(self,x,a,b)
    Computes the linear velocity field in 1 direction (scalar valued function) 
    such that :math:`u(x) = ax + b` for the given direction :math:`x` and 
    coefficients :math:`a` and :math:`b` computed by the method 
    :meth:`velocity_coefficients_1d`.

    :param x: coordinate of the direction in which the velocity varies
    :param a: slope of the component of the velocity
    :param b: constant of the component of the velocity

    :return: scalar valued function of the velocity field
    """
    return a*x + b
  
  def linear_velocity(self,x):
    """
    linear_velocity(self,x)
    computes the vector valued linear velocity function
    such that :math:`\\mathbf u(x) = \\mathbf a x + \\mathbf b`.
    Calls the method :meth:`linear_velocity_1d` for each direction.

    :param x: coordinates of the direction in which the velocity varies

    :return: vector valued velocity function
    """
    u = []
    for d in range(self.dim):
      u.append(self.linear_velocity_1d(x,self.a[d],self.b[d]))
    return u
  
  def evaluate_velocity_symbolic(self):
    """
    evaluate_velocity_symbolic(self)
    Evaluates the velocity field in symbolic form returning 
    a vector valued function of the form :math:`\\mathbf u (x,y,z) = \\mathbf a x + \\mathbf b y + \\mathbf c z + \\mathbf d`
    and shape ``(1,dim)``.
    Calls the method :meth:`linear_velocity` to evaluate the velocity field.
    """
    # create a numpy array of shape (1,dim) with the symbolic coordinates
    coor = np.array([[*self.sym_coor]], dtype='object')
    # rotate referential
    coor_R = self.rotate_referential(coor,self.O,self.L,ccw=False)
    # evaluate velocity with the rotated coordinate
    u = np.zeros(shape=(1,self.dim), dtype='object')
    u[0,:] = self.linear_velocity(coor_R[0,self.dir])
    # rotate the velocity field
    R   = self.rotation_matrix()
    u_R = self.rotate_vector(R,u,ccw=True)
    return u_R
  
  def evaluate_derivatives(self,u):
    """
    evaluate_derivatives(self,u)
    Evaluates the derivatives of the velocity function in symbolic form 
    returning a matrix of the shape ``(dim,dim)`` such that:
      
    .. math:: 
      \\nabla u_{ij} = \\frac{\\partial u_i}{\\partial x_j}

    
    :param u: vector valued function of the velocity field

    :return: **grad_u**: matrix of the derivatives of the velocity field shape ``(dim,dim)``
    """
    grad_u = np.zeros(shape=(self.dim,self.dim), dtype='object')
    for i in range(self.dim):
      for j in range(self.dim):
        grad_u[i,j] = u[i].diff(self.sym_coor[j])
    return grad_u

  def evaluate_velocity_and_derivatives_symbolic(self):
    """
    evaluate_velocity_and_derivatives_symbolic(self)
    Calls the methods 
    :meth:`evaluate_velocity_symbolic` 
    and 
    :meth:`evaluate_derivatives` 
    to evaluate the velocity field and its derivatives in symbolic form
    """
    u      = self.evaluate_velocity_symbolic()
    u[0,1] = self.evaluate_vertical_velocity(self.sym_coor[1],u=u)
    grad_u = self.evaluate_derivatives(u[0,:])
    return u,grad_u
  
  def evaluate_u_dot_n(self,u=None):
    """
    evaluate_u_dot_n(self,u=None)
    Evaluates the dot product of the velocity field with the normal vector in symbolic form 
    such that

    .. math:: 
      \\mathbf u \\cdot \\mathbf n = \\sum_{i=1}^{dim} u_i n_i
    
    :param u: **(Optional)**, vector valued function of the velocity field.
              If not provided, the symbolic velocity field is evaluated with :meth:`evaluate_velocity_symbolic`
    
    :return: **u_dot_n**: dot product of the velocity field with the normal vector
    """
    nmap = {1:'n_x',2:'n_x n_y',3:'n_x n_y n_z'}
    n = sp.symbols(nmap[self.dim])
    if u is None:
      u = self.evaluate_velocity_symbolic()
    
    u = sp.Matrix(u[0,:])
    n = sp.Matrix([n])
    u_dot_n = u.dot(n)
    return u_dot_n

  def evaluate_int_u_dot_n_faces(self,u=None):
    """
    evaluate_int_u_dot_n_faces(self,u=None)
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
    
    :param u: **(Optional)**, vector valued function of the velocity field.
              If not provided, the symbolic velocity field is evaluated with :meth:`evaluate_velocity_symbolic`

    :return: **int_u_dot_n**: dictionary with the integral of the dot product of the velocity field with the normal vector over the faces of the domain
    :rtype: dict
    """
    # evaluate symbolic u.n
    u_dot_n = self.evaluate_u_dot_n(u=u)
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

        # substitute the coordinate components by numerical values corresponding to the face
        _x = np.zeros(shape=(self.dim), dtype=np.float64)
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
        

        int_u_dot_n[face] = integral_u_dot_n
    return int_u_dot_n
  
  def evaluate_vertical_velocity_coefficients(self,u=None):
    """
    evaluate_vertical_velocity_coefficients(self,u=None)
    Evaluates the coefficients of the vertical component of the velocity function
    :math:`a` and :math:`b` such that :math:`u_y(y) = a y + b`
    To evaluate the coefficients, the following boundary conditions are used

    .. math::

      u_y(L_y) &= 0 \\\\
      u_y(O_y) &= \\frac{1}{S_{xz}} \\sum_f \\int_S \\mathbf u \\cdot \\mathbf n \\, dS \\\\

    where :math:`S_{xz}` is the surface of the bottom face of the domain, 
    :math:`O_y` and :math:`L_y` are the coordinates of the origin and max of the :math:`y` direction.

    :param u: **(Optional)**, vector valued function of the velocity field.
              If not provided, the symbolic velocity field is evaluated with :meth:`evaluate_velocity_symbolic`
    
    :return: **a**: slope of the vertical component of the velocity,
             **b**: constant of the vertical component of the velocity
    """
    if u is None:
      u = self.evaluate_velocity_symbolic()
    int_u_dot_n = self.evaluate_int_u_dot_n_faces(u=u)
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
    :meth:`evaluate_vertical_velocity(self,y,u=None)`
    Evaluates the vertical component of the velocity field

    :param y: coordinate of the vertical direction
    :param u: **(Optional)**, vector valued function of the velocity field.
              If not provided, the symbolic velocity field is evaluated with :meth:`evaluate_velocity_symbolic`

    :return: vertical component of the velocity field
    """
    if self.vertical_evaluated == False:
      self.evaluate_vertical_velocity_coefficients(u=u)
    uy = self.linear_velocity_1d(y,self.a[1],self.b[1])
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
    u = self.linear_velocity(coor_R[:,self.dir])
    # convert to numpy array of shape (npoints, dim)
    u = np.asarray(u, dtype=np.float64).T
    # rotate the velocity field
    R = self.rotation_matrix()
    u_R = self.rotate_vector(R,u,ccw=True)
    u_R[:,1] = self.evaluate_vertical_velocity(coor[:,1])
    return u_R
  
  def get_velocity_orientation(self,horizontal=True,normalize=True):
    """
    get_velocity_orientation(self,horizontal=True,normalize=False)
    Returns the orientation vector of the velocity field at the boundary

    :param bool horizontal: if True, only the horizontal components are returned (default: True)
    :param bool normalize:  if True, the vector is normalized (default: True)

    :return: **uL**: orientation of the velocity field at the boundary
    """
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