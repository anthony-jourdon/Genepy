#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: popctrl.py
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

class Markers:
  def __init__(self,layout:tuple[int,int,int]=(8,8,8)) -> None:
    self.layout = layout
    
class MarkersManagement(Markers):
  """
  .. py:class:: MarkersManagement(layout:tuple[int,int,int]=(8,8,8),popctrl_faces:tuple[int,...]=(0,1,4,5),popctrl_np_lower:int=8,popctrl_np_upper:int=128,popctrl_layout:tuple[int,int,int]=(2,2,2))

    Class to manage the lagrangian markers in `pTatin3d`_ model.

    :param layout: Number of markers per element in each direction. Default is (8,8,8).
    :type layout: tuple[int,int,int]
    :param popctrl_faces: Faces where markers are injected and not cleaned as long as they belong to an element connected to the face. Default is (0,1,4,5).
    :type popctrl_faces: tuple[int,...]
    :param popctrl_np_lower: Minimum number of markers per element. Default is 8.
    :type popctrl_np_lower: int
    :param popctrl_np_upper: Maximum number of markers per element. Default is 128.
    :type popctrl_np_upper: int
    :param popctrl_layout: Number of markers injected in each direction. Default is (2,2,2).
    :type popctrl_layout: tuple[int,int,int]
  """
  def __init__(self,layout:tuple[int,int,int]=(8,8,8),
               popctrl_faces:tuple[int,...]=(0,1,4,5),
               popctrl_np_lower:int=8,popctrl_np_upper:int=128,
               popctrl_layout:tuple[int,int,int]=(2,2,2)) -> None:
    Markers.__init__(self,layout)
    self.popctrl_faces    = popctrl_faces
    self.popctrl_np_lower = popctrl_np_lower
    self.popctrl_np_upper = popctrl_np_upper
    self.popctrl_layout   = popctrl_layout
    return
  
  def sprint_option(self):
    components = {"x":0,"y":1,"z":2}

    s = "###### Initial marker layout ######\n"
    for d in components:
      s += f"-lattice_layout_N{d} {self.layout[components[d]]} # markers per element in {d} direction\n"
    s += "###### Marker management ######\n"
    s += f"-mp_popctrl_np_lower {self.popctrl_np_lower} # min markers per element\n"
    s += f"-mp_popctrl_np_upper {self.popctrl_np_upper} # max markers per element\n"
    s += f"# marker injection in cells\n"
    for d in components:
      s += f"-mp_popctrl_n{d}p {self.popctrl_layout[components[d]]} # marker injection per element in {d} direction\n"
    s += f"# Boundary conditions for markers, perform injection and no cleaning of markers on marked faces\n"
    s += f"# Faces numbering:\n"
    s += f"#  0: east  = xmax = imax = Pxi\n"
    s += f"#  1: west  = xmin = imin = Nxi\n"
    s += f"#  2: north = ymax = jmax = Peta\n"
    s += f"#  3: south = ymin = jmin = Neta\n"
    s += f"#  4: front = zmax = kmax = Pzeta\n"
    s += f"#  5: back  = zmin = kmin = Nzeta\n"
    s += f"-model_GENE3D_bc_marker_nfaces {len(self.popctrl_faces)} # number of faces on which cleaning is ignored\n"
    s += f"-model_GENE3D_bc_marker_faces_list "
    for n in range(len(self.popctrl_faces)-1):
      s += f"{self.popctrl_faces[n]},"
    s += f"{self.popctrl_faces[len(self.popctrl_faces)-1]} # faces ignoring cleaning\n"
    return s
  
