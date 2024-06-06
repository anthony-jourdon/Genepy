#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: rotation.py
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

class Rotation:
  """
  .. py:class:: Rotation(dim,theta,axis=np.array([0,1,0]))

    Class to perform rotation of a referential in 2D or 3D.

    :param int dim: dimension of the rotation (can be 2 or 3)
    :param float theta: angle of rotation in radians
    :param np.ndarray axis: axis of rotation (default is y-axis)

    Example
    -------

    .. code-block:: python
      
        dim   = 3                                   # dimension
        theta = np.deg2rad(-90.0)                   # angle of rotation
        axis  = np.array([0,1,0], dtype=np.float64) # axis of rotation
        # Create class instance
        Rotation = genepy.Rotation(dim,theta,axis)

    Attributes
    ----------

    .. py:attribute:: dim
      :type: int
      :canonical: genepy.rotation.Rotation.dim

        Spatial dimension in which the rotation is performed
      
    .. py:attribute:: theta
      :type: float
      :canonical: genepy.rotation.Rotation.theta

        Angle of rotation in radians

    .. py:attribute:: axis
      :type: np.ndarray
      :canonical: genepy.rotation.Rotation.axis

        Axis of rotation, expected shape: ``(dim,)`` and dtype: ``np.float64``

    Methods
    -------
  """
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
    rotation_matrix_2d(self)
    computes the 2D rotation matrix given the angle :attr:`theta <genepy.rotation.Rotation.theta>`.
  
    :return: **R** rotation matrix in 2D of the shape ``(2,2)``
    """
    R = np.array([ [ np.cos(self.theta), -np.sin(self.theta) ],
                   [ np.sin(self.theta),  np.cos(self.theta) ] ], dtype=np.float64)
    return R
  
  def rotation_matrix_3d(self):
    """
    rotation_matrix_3d(self)
    computes the 3D rotation matrix given the angle :attr:`theta <genepy.rotation.Rotation.theta>` and the :attr:`axis <genepy.rotation.Rotation.axis>` of rotation.
  
    :return: **R** rotation matrix in 3D of the shape ``(3,3)``
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
    """
    rotation_matrix(self)
    Return the rotation matrix depending on the :attr:`spatial dimension <genepy.rotation.Rotation.dim>`.
    calls :meth:`rotation_matrix_2d() <genepy.rotation.Rotation.rotation_matrix_2d>` or :meth:`rotation_matrix_3d() <genepy.rotation.Rotation.rotation_matrix_3d>`.

    :return: **R** rotation matrix in 2D or 3D of the shape ``(2,2)`` or ``(3,3)``
    """
    if self.dim == 2:
      return self.rotation_matrix_2d()
    elif self.dim == 3:
      return self.rotation_matrix_3d()
    else:
      raise RuntimeError(f'Rotation can only be performed in 2D or 3D, found {self.dim}')
  
  def rotate_vector(self,R,u,ccw=True):
    """
    rotate_vector(self,R,u,ccw=True)
    Rotate vector(s) :math:`\\mathbf u` given the rotation matrix :math:`\\boldsymbol R`.

    .. warning:: 

      This is **not** a rotation of the vector field, but a rotation of the vectors themselves.
      To rotate the vector field, have a look at how it is done in
      :meth:`evaluate_velocity_numeric() <genepy.Velocity.evaluate_velocity_numeric>`.

    :param np.ndarray R: rotation matrix of the shape ``(dim,dim)``
    :param np.ndarray u: vector(s) to be rotated of the shape ``(npoints,dim)``
    :param bool ccw: rotate counter-clockwise (default is True)

    :return: **u_R** rotated vector(s) of the shape ``(npoints,dim)``
    """
    # u is expected to be in the form (npoints,dim)
    if ccw: # rotate couter-clockwise
      u_R = np.matmul(R,u.T).T
    else: # rotate clockwise
      u_R = np.matmul(R.T,u.T).T
    return u_R

  def rotate_referential(self,coor,O,L,ccw=True):
    """
    rotate_referential(self,coor,O,L,ccw=True)
    Rotate the referential of the coordinates :math:`\\mathbf{x}` given the rotation matrix :math:`\\boldsymbol R`.
    The referential is first translated to be centred on :math:`\\mathbf{0}`, then rotated and finally translated back to its original position.
    
    .. math:: 
      \\mathbf x_T &= \\mathbf x - \\frac{1}{2}(\\mathbf L + \\mathbf O) \\\\
      \\mathbf x_{TR} &= \\boldsymbol R \\mathbf x_T \\\\
      \\mathbf x_R &= \\mathbf x_{TR} + \\frac{1}{2}(\\mathbf L + \\mathbf O)

    :param np.ndarray coor: coordinates to be rotated of the shape ``(npoints,dim)``
    :param np.ndarray O: origin of the referential of the shape ``(dim,)``
    :param np.ndarray L: maximum coordinates of the referential of the shape ``(dim,)``
    :param bool ccw: rotate counter-clockwise (default is True)

    :return: **coorR** rotated coordinates of the shape ``(npoints,dim)``
    """
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
