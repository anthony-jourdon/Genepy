import numpy as np
import sympy as sp

class Rotation:
  def __init__(self,dim,theta,axis=np.array([0,1,0])) -> None:
    self.dim   = dim
    self.theta = theta
    self.axis  = axis
    return
  
  def __str__(self) -> str:
    attributes = vars(self)
    s = f'{self.__class__.__name__}:\n'
    for attribute in attributes:
      s += f'\t{attribute}:\t{attributes[attribute]}\n'
    return s
  
  def rotation_matrix_2d(self):
    """
    rotation_matrix_2d()
    computes the 2D rotation matrix given the angle theta.
  
    Parameters:
    -----------
    theta : angle of rotation
  
    Returns:
    --------
    R : rotation matrix
    """
    R = np.array([ [ np.cos(self.theta), -np.sin(self.theta) ],
                   [ np.sin(self.theta),  np.cos(self.theta) ] ], dtype=np.float64)
    return R
  
  def rotation_matrix_3d(self):
    """
    rotation_matrix_3d()
    computes the 3D rotation matrix given the angle theta and the axis of rotation.
  
    Parameters:
    -----------
    theta : angle of rotation
    axis  : axis of rotation
  
    Returns:
    --------
    R : rotation matrix
    """
    self.axis = self.axis/np.linalg.norm(self.axis)
    R = np.array([ [ np.cos(self.theta) + self.axis[0]**2*(1-np.cos(self.theta)), 
                     self.axis[0]*self.axis[1]*(1-np.cos(self.theta)) - self.axis[2]*np.sin(self.theta), 
                     self.axis[0]*self.axis[2]*(1-np.cos(self.theta)) + self.axis[1]*np.sin(self.theta) ],
                   [ self.axis[1]*self.axis[0]*(1-np.cos(self.theta)) + self.axis[2]*np.sin(self.theta), 
                     np.cos(self.theta) + self.axis[1]**2*(1-np.cos(self.theta)), 
                     self.axis[1]*self.axis[2]*(1-np.cos(self.theta)) - self.axis[0]*np.sin(self.theta) ],
                   [ self.axis[2]*self.axis[0]*(1-np.cos(self.theta)) - self.axis[1]*np.sin(self.theta), 
                     self.axis[2]*self.axis[1]*(1-np.cos(self.theta)) + self.axis[0]*np.sin(self.theta), 
                     np.cos(self.theta) + self.axis[2]**2*(1-np.cos(self.theta)) ] ], dtype=np.float64)
    return R
  
  def rotation_matrix(self):
    if self.dim == 2:
      return self.rotation_matrix_2d()
    elif self.dim == 3:
      return self.rotation_matrix_3d()
    else:
      raise RuntimeError(f'Rotation can only be performed in 2D or 3D, found {self.dim}')
  
  def rotate_vector(self,R,u,ccw=True):
    # u is expected to be in the form (npoints,dim)
    if ccw: # rotate couter-clockwise
      u_R = np.matmul(R,u.T).T
    else: # rotate clockwise
      u_R = np.matmul(R.T,u.T).T
    return u_R

  def rotate_referential(self,coor,O,L,ccw=True):
    # coor is expected to be in the form (npoints,dim)
    if coor.shape[1] != self.dim:
      raise RuntimeError(f'Coordinate dimension must be (npoint,{self.dim}), found {coor.shape}')
    # get the rotation matrix
    R = self.rotation_matrix()
    # translate referential to get centred on 0
    coorT = coor - 0.5*(L + O)
    # rotate
    coorTR = self.rotate_vector(R,coorT,ccw)
    # translate back
    coorR = coorTR + 0.5*(L + O)
    return coorR
  
def test2d():
  #O = np.array([-0.5, -0.5], dtype=np.float64)
  #L = np.array([0.5, 0.5], dtype=np.float64)
  O = np.array([0, 0], dtype=np.float64)
  L = np.array([1, 1], dtype=np.float64)
  x,z = sp.symbols('x z')
  coor = np.array([[x,z]], dtype='object')
  theta = np.deg2rad(-90.0)
  rot = Rotation(2,theta)
  print(rot)
  coor_R = rot.rotate_referential(coor,O,L,ccw=True)
  print(coor_R.shape)
  print('coor_R[0] =',coor_R[0,0])
  print('coor_R[1] =',coor_R[0,1])

def test3d():
  #O = np.array([-0.5, -0.5, -0.5], dtype=np.float64)
  #L = np.array([0.5, 0.5, 0.5], dtype=np.float64)
  O = np.array([0, 0, 0], dtype=np.float64)
  L = np.array([1, 1, 1], dtype=np.float64)
  x,y,z = sp.symbols('x y z')
  coor = np.array([[x,y,z]], dtype='object')
  theta = np.deg2rad(-90.0)
  axis = np.array([0,1,0], dtype=np.float64)
  rot = Rotation(3,theta,axis)
  print(rot)
  coor_R = rot.rotate_referential(coor,O,L,ccw=True)
  print(coor_R.shape)
  print('coor_R[0] =',coor_R[0,0])
  print('coor_R[1] =',coor_R[0,1])
  print('coor_R[2] =',coor_R[0,2])

if __name__ == "__main__":
  print("Running tests")
  test2d()
  test3d()