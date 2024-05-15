from . import domain
from . import rotation
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

class BoundaryConditions(domain.Domain,rotation.Rotation):
  def __init__(self,u_norm,u_angle,variation_dir:str,velocity_type:str,Domain,Rotation=None) -> None:
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

    self.velocity_coefficients()

  def __str__(self) -> str:
    attributes = vars(self)
    s  = f'{self.__class__.__name__}:\n'
    for attribute in attributes:
      if attribute == 'num_coor':
        s += f'\t{attribute}.shape:\t{attributes[attribute][0].shape}\n'
      else:
        s += f'\t{attribute}:\t{attributes[attribute]}\n'
    return s

  def boundary_vector(self):
    """
    boundary_vector(self)
    Computes the vector components based on the norm and angle with z axis (North-South) 
    Warning: these values are always positive, special attention is required if different signs are needed
    """
    uz = self.norm * np.cos(self.alpha)
    ux = np.sqrt(self.norm**2 - uz**2)
    return ux,uz
  
  def velocity_coefficients(self):
    """ 
    velocity_coefficients(self)
    computes velocity function coefficients such that
      u(x) = a*x + b
    with u = [ux, uz]
         a = [a0, a1]
         b = [b0, b1]
    """
    self.a = np.zeros(shape=(2), dtype=np.float64)
    self.b = np.zeros(shape=(2), dtype=np.float64)

    self.uO = np.zeros(shape=(2), dtype=np.float64)
    self.uL = np.zeros(shape=(2), dtype=np.float64)
    if self.type == "compression":
      self.uO[0], self.uO[1] = self.boundary_vector()
      self.uL = -self.uO
    elif self.type == "extension":
      self.uL[0], self.uL[1] = self.boundary_vector()
      self.uO = -self.uL
    else:
      raise RuntimeError(f'velocity_type can only be \"extension\" or \"compression\", found {self.type}')

    for d in range(2):
      self.a[d] = (self.uL[d] - self.uO[d]) / (self.L[self.dir] - self.O[self.dir])
      self.b[d] = -self.a[d]*self.L[self.dir] + self.uL[d]
    return

  def linear_velocity(self,x):
    ux = self.a[0]*x + self.b[0]
    uz = self.a[1]*x + self.b[1]
    return ux,uz
  
  def evaluate_velocity_symbolic(self):
    coor = np.array([[*self.sym_coor]], dtype='object')
    # rotate referential
    coor_R = self.rotate_referential(coor,self.O,self.L,ccw=False)
    # evaluate velocity with the rotated coordinate
    u = np.zeros(shape=(1,self.dim), dtype='object')
    u[0,0],u[0,self.dim-1] = self.linear_velocity(coor_R[0,self.dir])
    # rotate the velocity field
    R   = self.rotation_matrix()
    u_R = self.rotate_vector(R,u,ccw=True)
    return u_R
  
  def evaluate_derivatives(self,u):
    grad_u = np.zeros(shape=(self.dim,self.dim), dtype='object')
    for i in range(self.dim):
      for j in range(self.dim):
        grad_u[i,j] = u[i].diff(self.sym_coor[j])
    return grad_u

  def evaluate_velocity_and_derivatives_symbolic(self):
    u      = self.evaluate_velocity_symbolic()
    u[0,1] = self.evaluate_vertical_velocity(self.sym_coor[1],u=u)
    grad_u = self.evaluate_derivatives(u[0,:])
    return u,grad_u
  
  def evaluate_u_dot_n(self,u=None):
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

        # integrate u.n over the face simple integral if 2d, double integral if 3d
        if d == 0:
          integral_u_dot_n = sp.integrate(udn,(self.sym_coor[1],self.O[1],self.L[1]))
          if self.dim == 3: 
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
    
  def evaluate_vertical_velocity(self,y,u=None):
    if u is None:
      u = self.evaluate_velocity_symbolic()
    int_u_dot_n = self.evaluate_int_u_dot_n_faces(u=u)
    iudn = 0.0
    for face in int_u_dot_n:
      iudn += int_u_dot_n[face]
    bottom_surface = (self.L[0] - self.O[0])
    if self.dim == 3: bottom_surface *= (self.L[2] - self.O[2])
    uO = 0.0
    uL = -iudn / bottom_surface

    a = (uL - uO) / (self.L[1] - self.O[1])
    b = -a*self.L[1] + uL
    uy = a*y + b
    return uy
  
  def evaluate_velocity_numeric(self):
    coor = self.shape_coor()
    # rotate referential
    coor_R = self.rotate_referential(coor,self.O,self.L,ccw=False)
    # evaluate velocity with the rotated coordinate
    u = np.zeros(shape=(self.nv,self.dim), dtype=np.float64)
    u[:,0],u[:,self.dim-1] = self.linear_velocity(coor_R[:,self.dir])
    # rotate the velocity field
    R = self.rotation_matrix()
    u_R = self.rotate_vector(R,u,ccw=True)
    u_R[:,1] = self.evaluate_vertical_velocity(coor[:,1])
    return u_R
  
  def get_velocity_orientation(self,normalize=False):
    uL = np.zeros(shape=(self.dim), dtype=np.float64)
    uL[0],uL[self.dim-1] = self.uL
    R = self.rotation_matrix()
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
  
  def plot_velocity_vts(self,writer):
    writer.write_vts()
