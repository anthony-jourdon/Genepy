
import numpy as np
import sympy as sp

class Domain:
  """
  class Domain(dim,minCoor,maxCoor,size,coor=None)
  ------------------------------------------------

  Attributes:
  -----------
  dim:int          : dimension of the domain
  minCoor:np.array : minimum coordinates of the domain
  maxCoor:np.array : maximum coordinates of the domain
  size:np.array    : size of the domain (number of points in each direction)
  coor:np.array    : coordinates of the domain (numerical) if not provided, it will be created

  Methods:
  --------
  __init__(dim,minCoor,maxCoor,size,coor=None) : constructor
  __str__()                                    : string representation of the domain
  numerical_coordinates()                      : computes the numerical coordinates of the domain
  symbolic_coordinates()                       : computes the symbolic coordinates of the domain
  shape_coor()                                 : reshapes the coordinates from (nx,ny,nz) to (nv,dim
  """
  def __init__(self,dim:int,minCoor:np.ndarray,maxCoor:np.ndarray,size:np.ndarray,coor=None) -> None:
    self.dim   = dim
    
    if minCoor.shape[0] != self.dim:
      raise ValueError(f'{self.__class__.__name__}: minCoor must have {self.dim} elements, found {minCoor.shape[0]}')
    if maxCoor.shape[0] != self.dim:
      raise ValueError(f'{self.__class__.__name__}: maxCoor must have {self.dim} elements, found {maxCoor.shape[0]}')
    if size.shape[0] != self.dim:
      raise ValueError(f'{self.__class__.__name__}: size must have {self.dim} elements, found {size.shape[0]}')
    
    self.O  = minCoor # origin
    self.L  = maxCoor # length
    self.n  = size    # number of points in each direction
    self.nv = np.prod(self.n) # total number of points

    self.num_coor = coor # numerical coordinates
    self.sym_coor = None # symbolic coordinates
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
  
  def sprint_option(self,model_name:str):
    component = {0:'x', 1:'y', 2:'z'}
    s = "########### Bounding Box ###########\n"
    for d in range(self.dim):
      s += f"-{model_name}_O{component[d]} {self.O[d]:.5g} # min {component[d]} coord\n"
      s += f"-{model_name}_L{component[d]} {self.L[d]:.5g} # max {component[d]} coord\n"
    s += "########### Mesh ###########\n"
    for d in range(self.dim):
      s += f"-m{component[d]} {self.n[d]-1} # number of elements in {component[d]} direction\n"
    return s
  
  def sprint_option_list(self,model_name:str) -> list:
    component = {0:'x', 1:'y', 2:'z'}
    s = []
    s.append("########### Bounding Box ###########")
    for d in range(self.dim):
      s.append(f"-{model_name}_O{component[d]} {self.O[d]:.5g} # min {component[d]} coord")
      s.append(f"-{model_name}_L{component[d]} {self.L[d]:.5g} # max {component[d]} coord")
    s.append("########### Mesh ###########")
    for d in range(self.dim):
      s.append(f"-m{component[d]} {self.n[d]-1} # number of elements in {component[d]} direction")
    return s

  def numerical_coordinates(self):
    """
    numerical_coordinates(self)
    ---------------------------
    Computes the numerical coordinates of the domain. 
    Works in 1D, 2D and 3D.
    Attach the coordinates to the object as a tuple (X,Y,Z) with each 
    direction being of the shape X.shape = (nx,ny,nz).
    Tuples are immutable, so to modify the coordinates, convert them to a list first and restore them as tuple once done.
    """
    x1d = []
    for d in range(self.dim):
      x1d.append(np.linspace(self.O[d],self.L[d],self.n[d],dtype=np.float64))
    self.num_coor = tuple(np.meshgrid(*x1d,indexing='ij'))
    return
  
  def symbolic_coordinates(self):
    """
    symbolic_coordinates(self)
    -------------------------
    Computes the symbolic coordinates of the domain.
    Works in 1D, 2D and 3D.
    Attach the coordinates to the object as a tuple (x,y,z)
    """
    variables = {1:'x', 2:'x y', 3:'x y z'}
    self.sym_coor = sp.symbols(variables[self.dim])
    return
  
  def shape_coor(self):
    """
    shape_coor(self)
    ----------------
    Reshapes the coordinates from (nx,ny,nz) to (nv,dim)
    """
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