from . import domain
import numpy as np
import matplotlib.pyplot as plt

class BoundaryConditions(domain.Domain):
  def __init__(self,u_norm,u_angle,variation_dir:int,velocity_type:str,Domain) -> None:
    self.norm   = u_norm
    self.alpha  = u_angle
    self.dir    = variation_dir
    self.type   = velocity_type
    super(BoundaryConditions, self).__init__(Domain.O,Domain.L,Domain.n,Domain.theta)

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
  
  def rotate_vector(self,ux,uz,ccw=True):
    # rotate the vector field
    if ccw:
      ux_R = np.cos(self.theta) * ux - np.sin(self.theta) * uz
      uz_R = np.sin(self.theta) * ux + np.cos(self.theta) * uz
    else:
      ux_R =  np.cos(self.theta) * ux + np.sin(self.theta) * uz
      uz_R = -np.sin(self.theta) * ux + np.cos(self.theta) * uz
    return ux_R,uz_R

  def evaluate_velocity(self,x,z,ccw=True):
    xr,zr = self.rotate_referential(x,z)
    if self.dir == 0:
      coor = xr
    elif self.dir == 1:
      coor = zr
    # evaluate velocity with the rotated coordinate
    ux,uz   = self.linear_velocity(coor)
    uxr,uzr = self.rotate_vector(ux,uz,ccw=ccw)
    return uxr,uzr
  
  def evaluate_boundary_velocity(self):
    uO = self.rotate_vector(self.uO[0],self.uO[1],ccw=True)
    uL = self.rotate_vector(self.uL[0],self.uL[1],ccw=True)
    return uO,uL

  def plot_velocity(self,ccw=True):
    uxr,uzr = self.evaluate_velocity(self.num_coor[0],self.num_coor[1],ccw=ccw)
    u_norm = np.sqrt(uxr**2 + uzr**2)

    fig,ax = plt.subplots()
    ax.contourf(self.num_coor[0],self.num_coor[1],u_norm,100,cmap='Blues')
    ax.quiver(self.num_coor[0],self.num_coor[1],uxr,uzr)
    ax.axis('equal')
    plt.show()
    return
  
  def symbolic_velocity(self,ccw=True):
    # Evaluate velocity using symbols for x and z
    uxr,uzr = self.evaluate_velocity(self.sym_coor[0],self.sym_coor[1],ccw=ccw)
    print('###### velocity function ######')
    print('ux =',uxr)
    print('uz =',uzr)
    # Evaluate velocity at coords (Ox,Oz) and (Lx,Lz) 
    uO,uL = self.evaluate_boundary_velocity()
    print('###### velocity at (Ox,Oz) ######')
    print('uO =',uO)
    print('###### velocity at (Lx,Lz) ######')
    print('uL =',uL)
    return uxr,uzr
  
  def symbolic_derivatives(self,ux,uz):
    duxdx = ux.diff(self.sym_coor[0])
    duxdz = ux.diff(self.sym_coor[1])
    duzdx = uz.diff(self.sym_coor[0])
    duzdz = uz.diff(self.sym_coor[1])
    print('###### velocity derivatives ######')
    print('dux/dx =',duxdx)
    print('dux/dz =',duxdz)
    print('duz/dx =',duzdx)
    print('duz/dz =',duzdz)
    return

  def evaluate_velocity_and_derivatives(self,ccw=True):
    ux,uz = self.symbolic_velocity(ccw=ccw)
    self.symbolic_derivatives(ux,uz)
    return
