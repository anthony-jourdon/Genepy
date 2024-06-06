#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  bc-pre-processing
#  filename: pswarm.py
#
#  This file is part of bc-pre-processing.
#
#  bc-pre-processing is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, either 
#  version 3 of the License, or any later version.
#
#  bc-pre-processing is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with bc-pre-processing. 
#  If not, see <https://www.gnu.org/licenses/>.
#====================================================================================================

import numpy as np
from genepy import Markers

class Pswarm(Markers):
  """
  .. py:class:: Pswarm(layout:tuple[int,int,int]=(2,2,2),pressure:bool=False,temperature:bool=False,model_name:str="model_GENE3D")

    Class to define passive swarms of markers in `pTatin3d`_.
    It is the parent class of the different types of passive marker layout:

    - :py:class:`PswarmFillDomain`
    - :py:class:`PswarmFillDomainWithinBoundingBox`
    - :py:class:`PswarmFillBox`
    - :py:class:`PswarmFromUserList`

    :param layout: Number of markers per element in each direction. Default is (2,2,2).
    :type layout: tuple[int,int,int]
    :param pressure: Activate pressure tracking. Default is False.
    :type pressure: bool
    :param temperature: Activate temperature tracking. Default is False.
    :type temperature: bool
    :param model_name: Name of the model. Default is "model_GENE3D".
    :type model_name: str
  """
  def __init__(self, layout:tuple[int, int, int]=(2,2,2), 
               pressure:bool=False,
               temperature:bool=False,
               model_name:str="model_GENE3D") -> None:
    Markers.__init__(self,layout)
    self.prefix       = "passive_pswarm"
    self.coord_layout = -1
    self.model_name   = model_name
    self.pressure     = pressure
    self.temperature  = temperature

  def sprint_option(self) -> str:
    s  = "########### Passive tracers ###########\n"
    s += f"-{self.model_name}_{self.prefix}_apply # activate passive markers\n"
    s += f"-{self.model_name}_{self.prefix}_region_index 1\n"
    if self.pressure:
      s += f"-{self.model_name}_{self.prefix}_pressure # activate pressure tracking\n"
    if self.temperature:
      s += f"-{self.model_name}_{self.prefix}_temperature # activate temperature tracking\n"
    s += "###### Coordinate layout ######\n"
    s += "# 0: FillDM,  1: FillDMWithinBoundingBox\n"
    s += "# 2: FillBox, 3: FromUserList\n"
    s += f"-{self.model_name}_{self.prefix}_coord_layout {self.coord_layout} # type of marker filling\n"
    return s

  def __str__(self):
    s = "Passive swarm:\n"
    s += f"\tLayout: {self.layout}\n"
    s += f"\tTracking:\n"
    if self.pressure:
      s += f"\t\tPressure\n"
    if self.temperature:
      s += f"\t\tTemperature\n"
    return s
    
class PswarmFillDomain(Pswarm):
  """
  .. py:class:: PswarmFillDomain(layout:tuple[int,int,int]=(2,2,2),pressure:bool=False,temperature:bool=False,model_name:str="model_GENE3D")

    Class to define passive swarms of markers in `pTatin3d`_.
    Generates options for a passive swarm that fills the entire physical domain.

    :param layout: Number of markers per element in each direction. Default is (2,2,2).
    :type layout: tuple[int,int,int]
    :param pressure: Activate pressure tracking. Default is False.
    :type pressure: bool
    :param temperature: Activate temperature tracking. Default is False.
    :type temperature: bool
    :param model_name: Name of the model. Default is "model_GENE3D".
    :type model_name: str
  """
  def __init__(self, layout:tuple[int, int, int]=(2,2,2), 
               pressure:bool=False,
               temperature:bool=False,
               model_name:str="model_GENE3D") -> None:
    Pswarm.__init__(self,layout,pressure,temperature,model_name)
    self.coord_layout = 0
  
  def sprint_option(self) -> str:
    s = Pswarm.sprint_option(self)
    return s
  
  def __str__(self):
    s = Pswarm.__str__(self)
    s += "\tCoordinate layout: Fill Domain\n"
    return s
  
class PswarmFillDomainWithinBoundingBox(Pswarm):
  """
  .. py:class:: PswamFillDomainWithinBoundingBox(minCoord:tuple[float,float,float],maxCoord:tuple[float,float,float],layout:tuple[int,int,int]=(2,2,2),pressure:bool=False,temperature:bool=False,model_name:str="model_GENE3D")

    Class to define passive swarms of markers in `pTatin3d`_.
    Generates options for a passive swarm that fills the physical domain within a given bounding box.

    :param minCoord: Minimum coordinates of the bounding box.
    :type minCoord: tuple[float,float,float]
    :param maxCoord: Maximum coordinates of the bounding box.
    :type maxCoord: tuple[float,float,float]
    :param layout: Number of markers per element in each direction. Default is (2,2,2).
    :type layout: tuple[int,int,int]
    :param pressure: Activate pressure tracking. Default is False.
    :type pressure: bool
    :param temperature: Activate temperature tracking. Default is False.
    :type temperature: bool
    :param model_name: Name of the model. Default is "model_GENE3D".
    :type model_name: str
  """
  def __init__(self, minCoord:tuple[float,float,float],
               maxCoord:tuple[float,float,float], 
               layout:tuple[int, int, int]=(2,2,2), 
               pressure:bool=False,
               temperature:bool=False,
               model_name:str="model_GENE3D") -> None:
    Pswarm.__init__(self,layout,pressure,temperature,model_name)
    self.coord_layout = 1
    self.O = minCoord
    self.L = maxCoord
  
  def sprint_option(self) -> str:
    s  = Pswarm.sprint_option(self)
    s += f"-{self.model_name}_{self.prefix}_nx {self.layout[0]},{self.layout[1]},{self.layout[2]} # markers per element in each direction\n"
    s += f"-{self.model_name}_{self.prefix}_lattice_min {self.O[0]:g},{self.O[1]:g},{self.O[2]:g} # min coords of the bounding box\n"
    s += f"-{self.model_name}_{self.prefix}_lattice_max {self.L[0]:g},{self.L[1]:g},{self.L[2]:g} # max coords of the bounding box\n"
    return 
  
  def __str__(self):
    s = Pswarm.__str__(self)
    s += "\tCoordinate layout: Fill Domain Within Bounding Box\n"
    s += "\tBounding Box:\n"
    s += f"\t\tMin: {self.O}\n"
    s += f"\t\tMax: {self.L}\n"
    return s

class PswarmFillBox(Pswarm):
  """
  .. py:class:: PswarmFillBox(minCoord:tuple[float,float,float],maxCoord:tuple[float,float,float],layout:tuple[int,int,int]=(2,2,2),pressure:bool=False,temperature:bool=False,model_name:str="model_GENE3D")

    Class to define passive swarms of markers in `pTatin3d`_.
    Generates options for a passive swarm that fills a given bounding box with a specified marker resolution.

    .. warning:: 
      Opposed to the other classes, the layout is not the number of markers 
      per element in each direction but the total number of markers in each direction.

    :param minCoord: Minimum coordinates of the bounding box.
    :type minCoord: tuple[float,float,float]
    :param maxCoord: Maximum coordinates of the bounding box.
    :type maxCoord: tuple[float,float,float]
    :param layout: Number of markers in each direction. Default is (2,2,2).
    :type layout: tuple[int,int,int]
    :param pressure: Activate pressure tracking. Default is False.
    :type pressure: bool
    :param temperature: Activate temperature tracking. Default is False.
    :type temperature: bool
    :param model_name: Name of the model. Default is "model_GENE3D".
    :type model_name: str
  """
  def __init__(self, minCoord:tuple[float,float,float],
               maxCoord:tuple[float,float,float], 
               layout:tuple[int, int, int] = (2, 2, 2), 
               pressure:bool = False, 
               temperature:bool = False, 
               model_name:str = "model_GENE3D") -> None:
    Pswarm.__init__(self,layout,pressure,temperature,model_name)
    self.coord_layout = 2
    self.O = minCoord
    self.L = maxCoord
  
  def sprint_option(self) -> str:
    s = Pswarm.sprint_option(self)
    s += f"-{self.model_name}_{self.prefix}_box_nx {self.layout[0]},{self.layout[1]},{self.layout[2]} # number of markers (total, not per element) in each direction\n"
    s += f"-{self.model_name}_{self.prefix}_box_min {self.O[0]:g},{self.O[1]:g},{self.O[2]:g} # min coords of the bounding box\n"
    s += f"-{self.model_name}_{self.prefix}_box_max {self.L[0]:g},{self.L[1]:g},{self.L[2]:g} # max coords of the bounding box\n"
    return s
  
  def __str__(self):
    s = Pswarm.__str__(self)
    s += "\tCoordinate layout: Fill Box\n"
    s += "\tBox:\n"
    s += f"\t\tMin: {self.O}\n"
    s += f"\t\tMax: {self.L}\n"
    return s

class PswarmFromUserList(Pswarm):
  """
  .. py:class:: PswarmFromUserList(markers_coor:np.ndarray,layout:tuple[int,int,int]=(2,2,2),pressure:bool=False,temperature:bool=False,model_name:str="model_GENE3D")

    Class to define passive swarms of markers in `pTatin3d`_.
    Generates options for a list of passive tracers. 
    
    .. note:: 
      This layout should be used only when a small number of tracers are required
      and that their initial coordinates are known.
    
    :param markers_coor: Coordinates of the markers.
    :type markers_coor: np.ndarray
    :param pressure: Activate pressure tracking. Default is False.
    :type pressure: bool
    :param temperature: Activate temperature tracking. Default is False.
    :type temperature: bool
    :param model_name: Name of the model. Default is "model_GENE3D".
    :type model_name: str
  """
  def __init__(self, markers_coor:np.ndarray, 
               pressure: bool = False, 
               temperature: bool = False, 
               model_name: str = "model_GENE3D") -> None:
    layout: tuple[int, int, int] = (2, 2, 2) # not used for this type of layout
    Pswarm.__init__(self,layout,pressure,temperature,model_name)
    # coor should be of the shape (npoints,3)
    self.coord_layout = 3
    self.mcoor = markers_coor
  
  def sprint_option(self) -> str:
    s = Pswarm.sprint_option(self)
    npoints = self.mcoor.shape[0]
    component = {"x":0,"y":1,"z":2}
    s += f"-{self.model_name}_{self.prefix}_coor_n {npoints} # number of markers\n"
    for d in component:
      s += f"-{self.model_name}_{self.prefix}_coor_{d} "
      for i in range(npoints-1):
        s += f"{self.mcoor[i,component[d]]:g},"
      s += f"{self.mcoor[npoints-1,component[d]]:g} # {d} coord of each marker\n"
    return 

  def __str__(self):
    s = Pswarm.__str__(self)
    s += "\tCoordinate layout: From User List\n"
    s += "\tMarkers coordinates:\n"
    for i in range(self.mcoor.shape[0]):
      s += f"\t\t{self.mcoor[i,:]}\n"
    return s
