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