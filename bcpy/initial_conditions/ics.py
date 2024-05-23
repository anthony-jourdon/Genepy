from bcpy import Domain

class InitialConditions:
  """
  class InitialConditions
  -----------------------
  Abstract class for initial conditions

  Parameters:
  -----------
  Domain     : Domain object
  velocity   : velocity function (str or sympy expression)
  model_name : name of the model (default: model_GENE3D)
  **kwargs   : keyword arguments
               currently available keyword arguments are:
                - mesh_refinement : MeshRefinement object
                - initial_strain  : Gaussian object 
  
  Methods:
  --------
  sprint_option : returns the string to be written in the options file
  """
  def __init__(self, Domain:Domain, velocity, model_name:str="model_GENE3D", **kwargs) -> None:
    self.model_name      = model_name
    self.Domain          = Domain
    self.u               = velocity
    self.kwargs          = kwargs
    self.possible_kwargs = ["mesh_refinement","initial_strain"]

    for kw in self.kwargs:
      if kw not in self.possible_kwargs:
        serr = f"Keyword argument {kw} not recognized. Possible values are {self.possible_kwargs}"
        raise ValueError(serr)
    return
  
  def sprint_option(self) -> str:
    s = self.Domain.sprint_option(self.model_name)
    s += self.sprint_initial_velocity()
    for kw in self.possible_kwargs:
      if kw in self.kwargs:
        s += self.kwargs[kw].sprint_option(self.model_name)
    return s
  
  def sprint_initial_velocity(self) -> str:
    components = {0:"x", 1:"y", 2:"z"}
    prefix = "ic_velocity"
    s  = "########### Initial velocity field ###########\n"
    s += f"-{self.model_name}_{prefix}_ndir 3\n"
    s += f"-{self.model_name}_{prefix}_dir 0,1,2\n"
    for dim in components:
      u_string = str(self.u[0,dim])
      u_split  = u_string.split()
      nmembers = len(u_split)
      u_nospace = ""
      for j in range(nmembers):
        u_nospace += u_split[j]
      s += f"-{self.model_name}_{prefix}_expression_{dim} {u_nospace}\n"
    return s