from bcpy.initial_conditions import domain
from bcpy import rotation
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

class Velocity(domain.Domain,rotation.Rotation):
  """
  class Velocity(domain.Domain,rotation.Rotation)
  ---------------------------------------------------------
  Class to evaluate symbolic and numeric linear velocity function and its derivatives in space.

  Attributes:
  -----------
  norm   : velocity norm
  alpha  : angle of the velocity field with the z axis
  type   : type of velocity field (extension or compression)
  dir    : direction in which the velocity varies
  uO     : velocity at the component origin
  uL     : velocity at the component max value
  a      : slope of the of the velocity function
  b      : constant term of the velocity function
  vertical_evaluated : flag to check if the vertical velocity has been evaluated

  Methods:
  --------
  __init__(Domain,u_norm,variation_dir:str,velocity_type:str,u_angle=0.0,Rotation=None) : constructor
  __str__() : returns a string with the class name and attributes of the class
  report_symbolic_functions(u,grad_u,uL) : returns a string with the symbolic velocity function, its derivatives and the boundary velocity orientation
  boundary_vector() : computes the vector components
  velocity_coefficients_1d(uO,uL,O,L) : computes velocity function coefficients such that u(x) = a*x + b
  velocity_coefficients() : computes velocity function coefficients such that u(x) = a*x + b for all dimensions
  linear_velocity_1d(x,a,b) : computes the linear velocity field in 1 direction (scalar valued function) such that u(x) = a*x + b
  linear_velocity(x) : computes the linear velocity field (vector valued function) such that u(x) = a*x + b
  evaluate_velocity_symbolic() : evaluates the velocity field in symbolic form returning a vector valued function
  evaluate_derivatives(u) : evaluates the derivatives of the velocity function in symbolic form returning a matrix of the shape (dim,dim)
  evaluate_velocity_and_derivatives_symbolic() : evaluates the velocity field and its derivatives in symbolic form returning a vector valued function and a matrix of the shape (dim,dim)
  evaluate_u_dot_n(u=None) : evaluates the dot product of the velocity field with the normal vector in symbolic form
  evaluate_int_u_dot_n_faces(u=None) : evaluates the integral of the dot product of the velocity field with the normal vector over the faces
  evaluate_vertical_velocity_coefficients(u=None) : evaluates the coefficients of the vertical velocity component
  evaluate_vertical_velocity(y,u=None) : evaluates the vertical velocity component
  evaluate_velocity_numeric() : evaluates the velocity field in numeric form returning a vector valued function
  get_velocity_orientation(horizontal=True,normalize=False) : returns the orientation (vector) of the velocity field at the boundary
  plot_velocity_matplotlib() : plots the velocity field using matplotlib
  """
  def __init__(self,Domain:domain.Domain,u_norm,variation_dir:str,velocity_type:str,u_angle=0.0,Rotation:rotation.Rotation=None) -> None:
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
    self.prefix = "bc"
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
    ------------------------------------------
    Returns a string with the symbolic velocity function, its derivatives and the boundary velocity orientation.
    Can be used with the print() function.

    Parameters:
    -----------
    u      : vector valued function of the velocity field
    grad_u : matrix of the derivatives of the velocity field shape (dim,dim)
    uL     : orientation of the velocity field at the boundary

    Returns:
    --------
    s : string with the symbolic velocity function, its derivatives and the boundary velocity orientation. 
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
    ---------------------
    Computes the vector components
    1D & 2D: only the norm is used
    3D:      based on the norm and angle with z axis (North-South) 
    Warning: these values are always positive, special attention is required if different signs are needed

    Returns:
    --------
    u : boundary velocity, shape (dim,)
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
    velocity_coefficients_1d(self)
    ------------------------------
    computes velocity function coefficients such that
      u(x) = a*x + b
    for a given dimension

    Parameters:
    -----------
    uO : velocity at the component origin
    uL : velocity at the component max value
    O  : origin of the component
    L  : max value of the component

    Returns:
    --------
    a : slope of the component of the velocity
    b : constant of the component of the velocity
    """
    a = (uL - uO) / (L - O)
    b = -a*L + uL
    return a,b
  
  def velocity_coefficients(self):
    """ 
    velocity_coefficients(self)
    ---------------------------
    computes velocity function coefficients such that
      u(x) = a*x + b
    for all dimensions
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
    linear_velocity_1d(self)
    ------------------------
    computes the linear velocity field in 1 direction (scalar valued function) such that
      u(x) = a*x + b

    Parameters:
    -----------
    x : coordinates of the direction in which the velocity varies
    a : slope of the required component of the velocity
    b : constant of the required component of the velocity

    Returns:
    --------
    u : scalar velocity function
    """
    return a*x + b
  
  def linear_velocity(self,x):
    """
    linear_velocity(self)
    ---------------------
    computes the linear velocity field (vector valued function) such that
      u(x) = a*x + b

    Parameters:
    -----------
    x : coordinates of the direction in which the velocity varies

    Returns:
    --------
    u : vector velocity function
    """
    u = []
    for d in range(self.dim):
      u.append(self.linear_velocity_1d(x,self.a[d],self.b[d]))
    return u
  
  def evaluate_velocity_symbolic(self):
    """
    evaluate_velocity_symbolic(self)
    --------------------------------
    Evaluates the velocity field in symbolic form returning a vector valued function of the form
      u (x,y,z) = a*x + b*y + c*z + d
    and shape (1,dim)
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
    ----------------------------
    Evaluates the derivatives of the velocity function in symbolic form 
    returning a matrix of the shape (dim,dim) such that:
      grad_u[i,j] = du_i/dx_j
    
    Parameters:
    -----------
    u : vector valued function of the velocity field

    Returns:
    --------
    grad_u : matrix of the derivatives of the velocity field shape (dim,dim)
    """
    grad_u = np.zeros(shape=(self.dim,self.dim), dtype='object')
    for i in range(self.dim):
      for j in range(self.dim):
        grad_u[i,j] = u[i].diff(self.sym_coor[j])
    return grad_u

  def evaluate_velocity_and_derivatives_symbolic(self):
    """
    evaluate_velocity_and_derivatives_symbolic(self)
    ------------------------------------------------
    Evaluates the velocity field and its derivatives in symbolic form 
    returning a vector valued function of the form
      u (x,y,z) = a*x + b*y + c*z + d
    and a matrix of the shape (dim,dim) such that:
      grad_u[i,j] = du_i/dx_j

    Returns:
    --------
    u : vector valued function of the velocity field
    grad_u : matrix of the derivatives of the velocity field shape (dim,dim)
    """
    u      = self.evaluate_velocity_symbolic()
    u[0,1] = self.evaluate_vertical_velocity(self.sym_coor[1],u=u)
    grad_u = self.evaluate_derivatives(u[0,:])
    return u,grad_u
  
  def evaluate_u_dot_n(self,u=None):
    """
    evaluate_u_dot_n(self,u=None)
    Evaluates the dot product of the velocity field with the normal vector in symbolic form such that
      u.n = n_x*u_x + n_y*u_y + n_z*u_z
    
    Parameters:
    -----------
    (Optional)
    u : vector valued function of the velocity field
    
    Returns:
    --------
    u_dot_n : dot product of the velocity field with the normal vector
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
    if self.vertical_evaluated == False:
      self.evaluate_vertical_velocity_coefficients(u=u)
    uy = self.linear_velocity_1d(y,self.a[1],self.b[1])
    return uy
  
  def evaluate_velocity_numeric(self):
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
  
  def get_velocity_orientation(self,horizontal=True,normalize=False):
    """
    get_velocity_orientation(self,horizontal=True,normalize=False)
    Returns the orientation (vector) of the velocity field at the boundary

    Parameters:
    -----------
    horizontal : if True, only the horizontal components are returned
    normalize  : if True, the vector is normalized

    Returns:
    --------
    uL : orientation of the velocity field at the boundary
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