#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: utils.py
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

def compute_angle_between_points_2d(point_0,point_1):
  """
  compute_angle_between_points_2d(point_0,point_1)
  Computes the angle between two points.

  Parameters:
  -----------
  point_0 : coordinates of the first point
  point_1 : coordinates of the second point

  Returns:
  --------
  angle : angle between the two points
  """

  dx = np.abs(point_0[0] - point_1[0])
  dy = np.abs(point_0[1] - point_1[1])
  angle = 0.0
  if dx > 1.0e-6:
    angle = np.arctan(dy/dx)
  return angle

def x_centre_from_angle(z,angle,domain_centre):
  """
  x_centre_from_angle(notch_centre,angle,domain_centre)
  Computes the x coordinate of the centre of the notch given 
   - the z coordinate of the notch,
   - the angle
   - the centre of the domain.

  Parameters:
  -----------
  z             : z coordinate of the notch centre
  angle         : angle between the notch and the domain centre
  domain_centre : 2d coordinates of the domain centre

  Returns:
  --------
  x : x coordinate of the notch centre
  """
  x = np.tan(angle) * (z - domain_centre[1]) + domain_centre[0]
  return x

def newton_raphson(f,df,x0,tol=1.0e-6,max_iter=100,scaling=1.0,report=False,**kwargs):
  """
  newton_raphson(f,df,x0,tol=1.0e-6,max_iter=100)
  Newton-Raphson method to find the root of a function.

  Parameters:
  -----------
  f        : function to find the root of
  df       : derivative of the function
  x0       : initial guess
  tol      : tolerance
  max_iter : maximum number of iterations

  Returns:
  --------
  x : root of the function
  """
  if report: print("Newton-Raphson method started")
  res = 1e32
  it  = 0
  x   = x0
  while res > tol:
    if report: print(f"Iteration: {it}")
    x = x - f(x,**kwargs)/df(x,**kwargs)
    res = np.abs(f(x,**kwargs)*scaling)
    if report: print(f"Residual: {res:g} - x: {x:g}")
    it += 1
    if it > max_iter:
      print("Newton-Raphson method did not converge")
      break
  return x