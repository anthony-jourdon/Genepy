
import numpy as np
import sympy as sp

class Domain:
  def __init__(self,dim:int,minCoor:np.ndarray,maxCoor:np.ndarray,size:np.ndarray) -> None:
    self.dim   = dim
    
    if minCoor.shape[0] != self.dim:
      raise ValueError(f'minCoor must have {self.dim} elements, found {minCoor.shape[0]}')
    if maxCoor.shape[0] != self.dim:
      raise ValueError(f'maxCoor must have {self.dim} elements, found {maxCoor.shape[0]}')
    if size.shape[0] != self.dim:
      raise ValueError(f'size must have {self.dim} elements, found {size.shape[0]}')
    
    self.O     = minCoor
    self.L     = maxCoor
    self.n     = size

    self.num_coor = None
    self.sym_coor = None

    self.numerical_coordinates()
    self.symbolic_coordinates()

  def __str__(self) -> str:
    return f'Domain: {self.dim}D\nOrigin: {self.O}\nLength: {self.L}\nSize: {self.n}'

  def numerical_coordinates_2d(self):
    X,Z = np.meshgrid(np.linspace(self.O[0],self.L[0],self.n[0], dtype=np.float64),
                      np.linspace(self.O[1],self.L[1],self.n[1], dtype=np.float64), indexing='ij')
    self.num_coor = (X,Z)
    return
  
  def numerical_coordinates_3d(self):
    X,Y,Z = np.meshgrid(np.linspace(self.O[0],self.L[0],self.n[0], dtype=np.float64),
                        np.linspace(self.O[1],self.L[1],self.n[1], dtype=np.float64),
                        np.linspace(self.O[2],self.L[2],self.n[2], dtype=np.float64), indexing='ij')
    self.num_coor = (X,Y,Z)
    return

  def numerical_coordinates(self):
    if self.dim == 2:
      self.numerical_coordinates_2d()
    elif self.dim == 3:
      self.numerical_coordinates_3d()
    return
  
  def symbolic_coordinates_2d(self):
    self.sym_coor = sp.symbols('x z')
    return
  
  def symbolic_coordinates_3d(self):
    self.sym_coor = sp.symbols('x y z')
    return
  
  def symbolic_coordinates(self):
    if self.dim == 2:
      self.symbolic_coordinates_2d()
    elif self.dim == 3:
      self.symbolic_coordinates_3d()
    return

def test2d():
  O = np.array([-0.5, -0.5], dtype=np.float64)
  L = np.array([0.5, 0.5], dtype=np.float64)
  n = np.array([4,4], dtype=np.int32)
  dom = Domain(2,O,L,n)
  print(dom)
  print('X:',dom.num_coor[0].shape,'\n',dom.num_coor[0])
  print('Z:',dom.num_coor[1].shape,'\n',dom.num_coor[1])
  print(dom.sym_coor)
  coor = np.zeros(shape=(dom.n[0]*dom.n[1],2), dtype=np.float64)
  coor[:,0] = dom.num_coor[0].reshape(n[1]*n[0],order='F')
  coor[:,1] = dom.num_coor[1].reshape(n[1]*n[0],order='F')
  print(coor)

def test3d():
  O = np.array([-0.5, -0.5, -0.5], dtype=np.float64)
  L = np.array([0.5, 0.5, 0.5], dtype=np.float64)
  n = np.array([4,4,4], dtype=np.int32)
  dom = Domain(3,O,L,n)
  print(dom)
  print('X:',dom.num_coor[0].shape,'\n',dom.num_coor[0])
  print('Y:',dom.num_coor[1].shape,'\n',dom.num_coor[1])
  print('Z:',dom.num_coor[2].shape,'\n',dom.num_coor[2])
  print(dom.sym_coor)

  coor = np.zeros(shape=(dom.n[0]*dom.n[1]*dom.n[2],3), dtype=np.float64)
  coor[:,0] = dom.num_coor[0].reshape(n[2]*n[1]*n[0],order='F')
  coor[:,1] = dom.num_coor[1].reshape(n[2]*n[1]*n[0],order='F')
  coor[:,2] = dom.num_coor[2].reshape(n[2]*n[1]*n[0],order='F')
  print(coor)


if __name__ == "__main__":
  print("Running tests")
  test2d()
  print("\n")
  test3d()