class BoundaryCondition:
  def __init__(self, model_name:str="model_GENE3D") -> None:
    self.model_name = model_name

class StokesBoundaryCondition(BoundaryCondition):
  def __init__(self, tag:int=0, model_name:str="model_GENE3D") -> None:
    BoundaryCondition.__init__(self,model_name)
    self.tag:int        = tag
    self.prefix:str     = "bc"
    self.bc_type:int    = 0
    self.bc_name:str    = "uninitialized"
  
  def sprint_option(self):
    s  = f"###### Boundary condition {self.bc_name} {self.tag} ######\n"
    s += f"###### {self.__class__.__name__} boundary conditions ######\n"
    s += f"-{self.model_name}_bc_sc_name_{self.tag} {self.bc_name} # name of the boundary (arbitrary)\n"
    s += f"-{self.model_name}_bc_sc_type_{self.tag} {self.bc_type} # type of BC\n"
    s += f"-{self.model_name}_bc_facet_mesh_file_{self.tag} path_to_file # path to file of the boundary mesh\n"
    return s

  def __str__(self):
    s = f"{self.__class__.__name__}:\n"
    s += f"\tSurface tag: {self.tag}\n"
    s += f"\tName:        {self.bc_name}\n"
    s += f"\tType:        {self.bc_type}\n"
    return s

class TemperatureBC(BoundaryCondition):
  def __init__(self, faces:list[str], values:list[float], model_name:str="model_GENE3D") -> None:
    BoundaryCondition.__init__(self,model_name)
    self.faces  = faces
    self.values = values
  
  def sprint_option(self):
    s = f"###### Temperature boundary conditions ######\n"
    for face,value in zip(self.faces,self.values):
      if face not in ['xmin','xmax','ymin','ymax','zmin','zmax']:
        err = f"Error: face {face} not recognized.\n"
        err += "Recognized faces are: xmin, xmax, ymin, ymax, zmin, zmax.\n"
        raise ValueError(err)
      s += f"-{self.model_name}_bc_energy_{face} {value} # Temperature BC on face {face}\n"
    return s
  
  def __str__(self):
    s = "Temperature boundary conditions:\n"
    for face,value in zip(self.faces,self.values):
      s += f"\tFace: {face}\n"
      s += f"\tValue: {value}\n"
    return s

class ModelBCs:
  def __init__(self, velocity_bcs:list[StokesBoundaryCondition], 
               energy_bc:TemperatureBC=None, 
               model_name:str="model_GENE3D") -> None:
    self.model_name = model_name
    self.u_bcs      = velocity_bcs
    self.energy_bc  = energy_bc
    # enforce model name for all boundary conditions
    for bc in self.u_bcs:
      bc.model_name = self.model_name
    if self.energy_bc is not None:
      self.energy_bc.model_name = self.model_name

  def add_velocity_bc(self, bc:StokesBoundaryCondition) -> None:
    self.u_bcs.append(bc)
  
  def sprint_option(self):
    s = "########### Boundary conditions ###########\n"
    nsubfaces = len(self.u_bcs)
    s += f"-{self.model_name}_bc_nsubfaces {nsubfaces} # number of subfaces\n"
    s += f"-{self.model_name}_bc_tag_list "
    for n in range(nsubfaces-1):
      s += f"{self.u_bcs[n].tag},"
    s += f"{self.u_bcs[nsubfaces-1].tag} # list of boundary tags, determine the order of application\n"
    for bc in self.u_bcs:
      s += bc.sprint_option()
    if self.energy_bc is not None:
      s += self.energy_bc.sprint_option()
    return s

  def __str__(self):
    s = f"{self.model_name} boundary conditions:\n"
    for bc in self.u_bcs:
      s += str(bc)
    if self.energy_bc is not None:
      s += str(self.energy_bc)
    return s