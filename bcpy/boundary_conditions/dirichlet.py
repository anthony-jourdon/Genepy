from bcpy import StokesBoundaryCondition

class Dirichlet(StokesBoundaryCondition):
  def __init__(self, tag:int, name:str, components, velocity, model_name:str="model_GENE3D") -> None:
    StokesBoundaryCondition.__init__(self,tag,model_name)
    self.prefix     = "bc_dirichlet"
    self.bc_type    = 7
    self.bc_name    = name
    self.u          = velocity
    self.components = components
  
  def sprint_option(self):
    s = StokesBoundaryCondition.sprint_option(self)
    for d in self.components:
      if   d == "x": dim = 0
      elif d == "y": dim = 1
      elif d == "z": dim = 2
      u_string = str(self.u[0,dim])
      u_split  = u_string.split()
      nmembers = len(u_split)
      u_nospace = ""
      for j in range(nmembers):
        u_nospace += u_split[j]
      s += f"-{self.model_name}_{self.prefix}_u{d}_{self.tag} {u_nospace}\n"
    return s
  
  def __str__(self):
    s = StokesBoundaryCondition.__str__(self)
    s += f"\tComponents: {self.components}\n"
    s += f"\tVelocity:\n"
    for d in self.components:
      if   d == "x": dim = 0
      elif d == "y": dim = 1
      elif d == "z": dim = 2
      s += f"\t\tu{d}: {self.u[0,dim]}\n"
    return s

class DirichletUdotN(StokesBoundaryCondition):
  def __init__(self, tag:int, name:str, model_name:str="model_GENE3D") -> None:
    StokesBoundaryCondition.__init__(self,tag,model_name)
    self.prefix     = "bc_dirichlet"
    self.bc_type    = 7
    self.bc_name    = name
    
  def sprint_option(self):
    s = StokesBoundaryCondition.sprint_option(self)
    s += f"-{self.model_name}_{self.prefix}_bot_u.n_{self.tag}\n"
    return s
  
  def __str__(self):
    s = StokesBoundaryCondition.__str__(self)
    s += f"\tVelocity: balancing inflow/outflow\n"
    s += f"\ton bottom face using int u.n dS\n"
    return s