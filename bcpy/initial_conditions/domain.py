import numpy as np
import sympy as sp

class Domain:
  """
  .. py:class:: Domain(dim, minCoor, maxCoor, size, coor=None)

    Class to build a physical domain. 
    Compatible with 1D, 2D and 3D domains.
    Coordinates (symbolic and numeric) are instantiated by the class constructor using provided coordinates or creating new if not provided.

    :param int dim: dimension of the domain (can be 1, 2 or 3)
    :param np.ndarray minCoor: minimum coordinates of the domain
    :param np.ndarray maxCoor: maximum coordinates of the domain
    :param np.ndarray size: size of the domain (number of points in each direction)
    :param np.ndarray coor: coordinates of the domain (numerical), optional

    Example
    -------

    .. code-block:: python

      dim      = 3
      minCoord = np.array([0,-1,-0.5], dtype=np.float64)
      maxCoord = np.array([1,0,0.5], dtype=np.float64)
      size     = np.array([9,5,17], dtype=np.int32)

      domain   = bcpy.Domain(dim,minCoord,maxCoord,size)

    Attributes
    ----------

    .. py:attribute:: dim
      :type: int
      :canonical: bcpy.initial_conditions.domain.Domain.dim

        Dimension of the domain
    
    .. py:attribute:: O
      :type: np.ndarray
      :canonical: bcpy.initial_conditions.domain.Domain.O

      Minimum coordinates of the domain, expected shape: ``(dim,)`` and dtype: ``np.float64``

    .. py:attribute:: L
      :type: np.ndarray
      :canonical: bcpy.initial_conditions.domain.Domain.L

      Maximum coordinates of the domain, expected shape: ``(dim,)`` and dtype: ``np.float64``
    
    .. py:attribute:: n
      :type: np.ndarray
      :canonical: bcpy.initial_conditions.domain.Domain.n

      Size of the domain (number of nodes in each direction), expected shape: ``(dim,)`` and dtype: ``np.int32``

    .. py:attribute:: nv
      :type: int
      :canonical: bcpy.initial_conditions.domain.Domain.nv

      Total number of nodes in the domain

    .. py:attribute:: num_coor
      :type: tuple
      :canonical: bcpy.initial_conditions.domain.Domain.num_coor

      Numerical coordinates of the domain created by :py:meth:`numerical_coordinates() <bcpy.Domain.numerical_coordinates>` . 
      Tuple ``(X)``, ``(X,Y)`` or ``(X,Y,Z)`` (depending on the number of dimensions)
      with each direction being of type ndarray and shape ``X.shape = (*n)``.

    .. py:attribute:: sym_coor
      :type: tuple
      :canonical: bcpy.initial_conditions.domain.Domain.sym_coor

      Symbolic coordinates of the domain. 
      Tuple ``('x')``, ``('x','y')`` or ``('x','y','z')``

    Methods
    -------
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
    """
    __str__(self)
    Returns a string representation of the domain.
    """
    attributes = vars(self)
    s = f'{self.__class__.__name__}:\n'
    for attribute in attributes:
      if attribute == 'num_coor':
        s += f'\t{attribute}.shape:\t{attributes[attribute][0].shape}\n'
      else:
        s += f'\t{attribute}:\t{attributes[attribute]}\n'
    return s
  
  def sprint_option(self,model_name:str):
    """
    sprint_option(self,model_name:str)
    Returns a string formatted for `pTatin3d`_ input file using `PETSc <https://petsc.org>`_ options format.
    
    :param str model_name: name of the model to include in the options

    :return: string with the options
    """
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
    Computes the numerical coordinates of the domain as a uniform grid. 
    Compatible with 1D, 2D and 3D.
    Attach the coordinates to the attribute 
    :attr:`num_coor <bcpy.initial_conditions.domain.Domain.num_coor>` 
    as a tuple: ``self.num_coor = (X,Y,Z)`` 
    with each direction being of the shape ``X.shape = (self.n[0],self.n[1],self.n[2])``.
    Tuples are immutable, so to modify the coordinates, convert them to a list first and restore them as tuple once done.
    """
    x1d = []
    for d in range(self.dim):
      x1d.append(np.linspace(self.O[d],self.L[d],self.n[d],dtype=np.float64))
    self.num_coor = tuple(np.meshgrid(*x1d,indexing='ij'))
    return
  
  def symbolic_coordinates(self):
    """
    Computes the symbolic coordinates of the domain.
    Works in 1D, 2D and 3D.
    Attach the coordinates to the attribute 
    :py:attr:`sym_coor <bcpy.initial_conditions.domain.Domain.sym_coor>` 
    as a tuple (``"x"``, ``"y"``, ``"z"``)
    """
    variables = {1:'x', 2:'x y', 3:'x y z'}
    self.sym_coor = sp.symbols(variables[self.dim])
    return
  
  def shape_coor(self):
    """
    shape_coor(self)
    Reshapes the coordinates from ``(n[0],n[1],n[2])`` to ``(nv,dim)``
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