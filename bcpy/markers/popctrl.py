class Markers:
  def __init__(self,layout:tuple[int,int,int]=(8,8,8)) -> None:
    self.layout = layout
    
class MarkersManagement(Markers):
  """
  class MarkersManagement(Markers)
  --------------------------------
  MarkersManagement class is used to generate the options to manage the lagrangian markers in pTatin3d model.

  Attributes:
  -----------
  layout: tuple[int,int,int]
    Number of markers per element in each direction.
    Default: (8,8,8)
  popctrl_faces: tuple[int,...]
    Faces where markers are injected and not cleaned as long as they belong to an element connected to the face.
    Default: (0,1,4,5) -> east, west, front, back
  popctrl_np_lower: int
    Minimum number of markers per element.
    Default: 8
  popctrl_np_upper: int
    Maximum number of markers per element.
    Default: 128
  popctrl_layout: tuple[int,int,int]
    Number of markers injected in each direction.
    Default: (2,2,2)
  
  Methods:
  --------
  sprint_option():
    Returns the string to be included in the input file.
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
  