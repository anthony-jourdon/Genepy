
import numpy as np
import sympy as sp

class Domain:
  def __init__(self,minCoor=np.zeros(shape=(2), dtype=np.float64),
                    maxCoor=np.ones(shape=(2), dtype=np.float64),
                    size=np.array([2, 2], dtype=np.int32),
                    referential_angle=0.0) -> None:
    self.O     = minCoor
    self.L     = maxCoor
    self.n     = size
    self.theta = referential_angle

    self.num_coor = None
    self.sym_coor = None

    self.numerical_coordinates()
    self.symbolic_coordinates()

  def numerical_coordinates(self):
    coor = list()
    for i in range(self.n.shape[0]):
      coor.append( np.linspace(self.O[i],self.L[i],self.n[i], dtype=np.float64) )
    X,Z = np.meshgrid(coor[0],coor[1])
    self.num_coor = (X,Z)
    return

  def symbolic_coordinates(self):
    self.sym_coor = sp.symbols('x z')
    return
  
  def rotate_referential(self,x,z,ccw=False):
    # translate
    xT = x - 0.5*(self.L[0] + self.O[0])
    zT = z - 0.5*(self.L[1] + self.O[1])
    # rotate
    if ccw:
      xTR = xT * np.cos(self.theta) - zT * np.sin(self.theta)
      zTR = zT * np.cos(self.theta) + xT * np.sin(self.theta)
    else:
      xTR = xT * np.cos(self.theta) + zT * np.sin(self.theta)
      zTR = zT * np.cos(self.theta) - xT * np.sin(self.theta)
    # translate back
    xR = xTR + 0.5*(self.L[0] + self.O[0])
    zR = zTR + 0.5*(self.L[1] + self.O[1])
    return xR,zR