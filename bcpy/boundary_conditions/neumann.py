from bcpy import StokesBoundaryCondition

class Neumann(StokesBoundaryCondition):
  def __init__(self, tag:int, name:str, dev_stress_expression=None, mesh_file:str="path_to_file", model_name:str="model_GENE3D") -> None:
    StokesBoundaryCondition.__init__(self,tag,mesh_file,model_name)
    self.prefix     = "bc_neumann"
    self.bc_type    = 1
    self.bc_name    = name
    self.dev_stress = dev_stress_expression
  
  def sprint_option(self):
    s = StokesBoundaryCondition.sprint_option(self)
    if self.dev_stress is not None:
      s += f"-{self.model_name}_{self.prefix}_dev_stress_{self.tag} {self.dev_stress} # deviatoric stress\n"
    return s
  
  def __str__(self):
    s = StokesBoundaryCondition.__str__(self)
    s += f"Pressure component is automatically calculated\n"
    s += f"\tDeviatoric stress: {self.dev_stress}\n"
    return s