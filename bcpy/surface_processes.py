class SPM:
  def __init__(self, model_name:str="model_GENE3D") -> None:
    self.model_name = model_name
    self.prefix = "spm"

class SPMDiffusion(SPM):
  """
  class SPMDiffusion
  ------------------
  Surface processes: Diffusion computed in pTatin3d as
    \\nabla \\cdot (k \\nabla h) = 0
  with k the diffusivity and h the height.

  Attributes:
  -----------
  dirichlet_faces: list[str]
    List of faces where Dirichlet boundary conditions are applied.
    Recognized faces are: \"xmin\", \"xmax\", \"zmin\", \"zmax\".
    At least one face must be passed.
  diffusivity: float
    Diffusivity value.
  model_name: str (optional, default="model_GENE3D")
    Model name.

  Methods:
  --------
  sprint_option():
    Returns the string to be included in the input file.
  __str__():
    Returns the string representation of the object.
  """
  def __init__(self, dirichlet_faces:list[str],
               diffusivity:float=1.0e-6, 
               model_name:str="model_GENE3D") -> None:
    SPM.__init__(self,model_name)
    for face in dirichlet_faces:
      if face not in ['xmin','xmax','zmin','zmax']:
        err  = f"Error: SPM diffusivity Dirichlet boundary conditions:\n"
        err += f"Face {face} not recognized.\n"
        err += "Recognized faces are: xmin, xmax, zmin, zmax.\n"
        raise ValueError(err)
    self.bcs         = dirichlet_faces
    self.diffusivity = diffusivity
  
  def sprint_option(self):
    s = "########### Surface processes ###########\n"
    s += f"-{self.model_name}_{self.prefix}_apply_surface_diffusion # activate diffusion at surface\n"
    s += f"-{self.model_name}_{self.prefix}_diffusivity {self.diffusivity:g}\n"
    s += "# spm diffusion boundary conditions minimum 1 must be passed\n"
    for face in self.bcs:
      s += f"-{self.model_name}_{self.prefix}_diffusion_dirichlet_{face}\n"
    return s
  
  def __str__(self):
    s = "Surface processes: Diffusion\n"
    s += f"\tDiffusivity: {self.diffusivity}\n"
    s += f"\tDirichlet faces: {self.bcs}\n"
    return s