
import numpy as np
import sympy as sp

class Domain:
  def __init__(self,dim:int,minCoor:np.ndarray,maxCoor:np.ndarray,size:np.ndarray,coor=None) -> None:
    self.dim   = dim
    
    if minCoor.shape[0] != self.dim:
      raise ValueError(f'{self.__class__.__name__}: minCoor must have {self.dim} elements, found {minCoor.shape[0]}')
    if maxCoor.shape[0] != self.dim:
      raise ValueError(f'{self.__class__.__name__}: maxCoor must have {self.dim} elements, found {maxCoor.shape[0]}')
    if size.shape[0] != self.dim:
      raise ValueError(f'{self.__class__.__name__}: size must have {self.dim} elements, found {size.shape[0]}')
    
    self.O  = minCoor
    self.L  = maxCoor
    self.n  = size
    self.nv = np.prod(self.n)

    self.num_coor = coor
    self.sym_coor = None
    if self.num_coor is None: 
      self.numerical_coordinates()
    self.symbolic_coordinates()

  def __str__(self) -> str:
    attributes = vars(self)
    s = f'{self.__class__.__name__}:\n'
    for attribute in attributes:
      if attribute == 'num_coor':
        s += f'\t{attribute}.shape:\t{attributes[attribute][0].shape}\n'
      else:
        s += f'\t{attribute}:\t{attributes[attribute]}\n'
    return s

  def numerical_coordinates(self):
    x1d = []
    for d in range(self.dim):
      x1d.append(np.linspace(self.O[d],self.L[d],self.n[d],dtype=np.float64))
    self.num_coor = tuple(np.meshgrid(*x1d,indexing='ij'))
    return
  
  def symbolic_coordinates(self):
    variables = {1:'x', 2:'x y', 3:'x y z'}
    self.sym_coor = sp.symbols(variables[self.dim])
    return
  
  def shape_coor(self):
    coor = np.zeros(shape=(self.nv,self.dim), dtype=np.float64)
    for d in range(self.dim):
      coor[:,d] = self.num_coor[d].reshape(self.nv,order='F')
    return coor

def test2d():
  O = np.array([-0.5, -0.5], dtype=np.float64)
  L = np.array([0.5, 0.5], dtype=np.float64)
  n = np.array([4,3], dtype=np.int32)
  dom = Domain(2,O,L,n)
  print(dom)
  
def test3d():
  O = np.array([-0.5, -0.5, -0.5], dtype=np.float64)
  L = np.array([0.5, 0.5, 0.5], dtype=np.float64)
  n = np.array([4,5,3], dtype=np.int32)
  dom = Domain(3,O,L,n)
  print(dom)

if __name__ == "__main__":
  print("Running tests")
  test2d()
  print("\n")
  test3d()