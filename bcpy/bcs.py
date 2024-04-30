from . import domain
from . import rotation
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

class BoundaryConditions(domain.Domain,rotation.Rotation):
  def __init__(self,u_norm,u_angle,variation_dir:int,velocity_type:str,Domain,Rotation) -> None:
    self.norm   = u_norm
    self.alpha  = u_angle
    self.dir    = variation_dir
    self.type   = velocity_type

    domain.Domain.__init__(self,Domain.dim,Domain.O,Domain.L,Domain.n)
    rotation.Rotation.__init__(self,Domain.dim,Rotation.theta,Rotation.axis)

    self.velocity_coefficients()

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
  
  def evaluate_velocity_symbolic(self,x,y,z):
    if self.dim == 2:
      coor = np.array([[x,z]], dtype='object')
    elif self.dim == 3:
      coor = np.array([[x,y,z]], dtype='object')
    # rotate referential
    coor_R = self.rotate_referential(coor,self.O,self.L,ccw=False)
    # evaluate velocity with the rotated coordinate
    ux,uz = self.linear_velocity(coor_R[0,self.dir])
    # rotate the velocity field
    if self.dim == 2:
      u = np.array([[ux,uz]], dtype='object')
    elif self.dim == 3:
      u = np.array([[ux,0,uz]], dtype='object')
    R = self.rotation_matrix()
    u_R = self.rotate_vector(R,u,ccw=True)
    return u_R
  
  def evaluate_derivatives_2d(self,ux,uz,x,z):
    duxdx = ux.diff(x); duxdz = ux.diff(z)
    duzdx = uz.diff(x); duzdz = uz.diff(z)
    grad_u = sp.Matrix([[duxdx,duxdz],
                        [duzdx,duzdz]])
    return grad_u

  def evaluate_derivatives_3d(self,ux,uy,uz,x,y,z):
    duxdx = ux.diff(x); duxdy = ux.diff(y); duxdz = ux.diff(z)
    duydx = uy.diff(x); duydy = uy.diff(y); duydz = uy.diff(z)
    duzdx = uz.diff(x); duzdy = uz.diff(y); duzdz = uz.diff(z)
    grad_u = sp.Matrix([[duxdx,duxdy,duxdz],
                        [duydx,duydy,duydz],
                        [duzdx,duzdy,duzdz]])
    return grad_u

  def evaluate_velocity_and_derivatives_symbolic(self):
    if self.dim == 2:
      u = self.evaluate_velocity_symbolic(self.sym_coor[0],None,self.sym_coor[1])
      ux,uz = u[0,0],u[0,1]
      grad_u = self.evaluate_derivatives_2d(ux,uz,self.sym_coor[0],self.sym_coor[1])
    elif self.dim == 3:
      u = self.evaluate_velocity_symbolic(self.sym_coor[0],self.sym_coor[1],self.sym_coor[2])
      ux,uy,uz = u[0,0],u[0,1],u[0,2]
      grad_u = self.evaluate_derivatives_3d(ux,uy,uz,self.sym_coor[0],self.sym_coor[1],self.sym_coor[2])
    return u,grad_u
  
  def evaluate_velocity_numeric(self):
    if self.dim == 2:
      coor = np.zeros(shape=(self.n[0]*self.n[1],2), dtype=np.float64)
      coor[:,0] = self.num_coor[0].reshape(self.n[1]*self.n[0],order='F')
      coor[:,1] = self.num_coor[1].reshape(self.n[1]*self.n[0],order='F')
    elif self.dim == 3:
      coor = np.zeros(shape=(self.n[0]*self.n[1]*self.n[2],3), dtype=np.float64)
      coor[:,0] = self.num_coor[0].reshape(self.n[2]*self.n[1]*self.n[0],order='F')
      coor[:,1] = self.num_coor[1].reshape(self.n[2]*self.n[1]*self.n[0],order='F')
      coor[:,2] = self.num_coor[2].reshape(self.n[2]*self.n[1]*self.n[0],order='F')
    # rotate referential
    coor_R = self.rotate_referential(coor,self.O,self.L,ccw=False)
    # evaluate velocity with the rotated coordinate
    ux,uz = self.linear_velocity(coor_R[:,self.dir])
    # rotate the velocity field
    if self.dim == 2:
      u = np.zeros(shape=(self.n[0]*self.n[1],2), dtype=np.float64)
      u[:,0] = ux
      u[:,1] = uz
    elif self.dim == 3:
      u = np.zeros(shape=(self.n[0]*self.n[1]*self.n[2],3), dtype=np.float64)
      u[:,0] = ux
      u[:,2] = uz
    R = self.rotation_matrix()
    u_R = self.rotate_vector(R,u,ccw=True)
    return u_R

  def plot_velocity_matplotlib(self):
    u = self.evaluate_velocity_numeric()
    if self.dim == 2:
      X = self.num_coor[0]
      Z = self.num_coor[1]
      ux = u[:,0].reshape(self.n[0],self.n[1])
      uz = u[:,1].reshape(self.n[0],self.n[1])
    elif self.dim == 3:
      X = self.num_coor[0]
      Z = self.num_coor[2]
      ux = u[:,0].reshape(self.n[0],self.n[1],self.n[2])
      uz = u[:,2].reshape(self.n[0],self.n[1],self.n[2])
    u_norm = np.sqrt(ux**2 + uz**2)

    fig,ax = plt.subplots()
    ax.contourf(X,Z,u_norm,100,cmap='Blues')
    ax.quiver(X,Z,ux,uz)
    ax.axis('equal')
    plt.show()
    return
  
  def plot_velocity_vts(self,writer):
    writer.write_vts()
