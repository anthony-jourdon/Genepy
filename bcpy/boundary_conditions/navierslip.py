from bcpy import StokesBoundaryCondition

class NavierSlip(StokesBoundaryCondition):
  def __init__(self, tag:int, name:str, grad_u, u_orientation, mesh_file:str="path_to_file", model_name:str="model_GENE3D") -> None:
    StokesBoundaryCondition.__init__(self,tag,mesh_file,model_name)
    self.prefix     = "bc_navier"
    self.bc_type    = 6
    self.bc_name    = name
    self.grad_u     = grad_u
    self.uL         = u_orientation
    
  def sprint_option(self):
    s = StokesBoundaryCondition.sprint_option(self)
    s += f"-{self.model_name}_{self.prefix}_penalty_{self.tag} 1.0e3 # penalty for Nitsche's method\n"
    s += f"-{self.model_name}_{self.prefix}_duxdx_{self.tag} {str(self.grad_u[0,0])}\n"
    s += f"-{self.model_name}_{self.prefix}_duxdz_{self.tag} {str(self.grad_u[0,2])}\n"
    s += f"-{self.model_name}_{self.prefix}_duzdx_{self.tag} {str(self.grad_u[2,0])}\n"
    s += f"-{self.model_name}_{self.prefix}_duzdz_{self.tag} {str(self.grad_u[2,2])}\n"
    s += f"-{self.model_name}_{self.prefix}_uL_{self.tag} {self.uL[0]},{self.uL[2]} # velocity orientation on boundary\n"
    s += f"-{self.model_name}_{self.prefix}_mathcal_H_{self.tag} 0,1,0,1,1,1\n"
    return s
  
  def __str__(self):
    s = StokesBoundaryCondition.__str__(self)
    s += f"\tVelocity gradient:\n"
    components = ['x','y','z']
    for i in range(3):
      for j in range(3):
        s += f"\t\tdu{components[i]}/d{components[j]}: {str(self.grad_u[i,j])}\n"
    s += f"\tOrientation of velocity: {self.uL[0]},{self.uL[2]}\n"
    return s