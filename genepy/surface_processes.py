#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: surface_processes.py
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

class SPM:
  def __init__(self, model_name:str="model_GENE3D") -> None:
    self.model_name = model_name
    self.prefix = "spm"

class SPMDiffusion(SPM):
  """
  .. py:class:: SPMDiffusion(dirichlet_faces:list[str], diffusivity:float=1.0e-6, model_name:str="model_GENE3D")

    Class to generate the options for surface processes diffusion in `pTatin3d`_.
    Surface processes are handled by solving the following diffusion equation:
    
    .. math:: 
      \\nabla \\cdot (k \\nabla h) = 0
    
    with :math:`k` the diffusivity and :math:`h` the height.

    :param dirichlet_faces: List of faces where Dirichlet boundary conditions are applied. 
                            Possible values are: ``"xmin"``, ``"xmax"``, ``"zmin"``, ``"zmax"``.
                            At least one is required to ensure uniqueness of the solution.
    :type dirichlet_faces: list[str]
    :param diffusivity: Diffusivity value in m\\ :sup:`2`.s\\ :sup:`-1`. Default is :math:`10^{-6}`.
    :type diffusivity: float
    :param model_name: Model name.
    :type model_name: str
    :raises ValueError: If an unrecognized face is passed.
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

class SPMDiffusionBaselvl(SPM):
  """
  .. py:class:: SPMDiffusionBaselvl(dirichlet_faces:list[str], diffusivity:list[float]=[1.0e-10,1.0e-6], baselevel:float=0.0, model_name:str="model_GENE3D")

    Class to generate the options for surface processes diffusion in `pTatin3d`_ with a base level.
    The equation solved is the same than in :py:class:`SPMDiffusion`.
    
    :param dirichlet_faces: List of faces where Dirichlet boundary conditions are applied. 
                            Possible values are: ``"xmin"``, ``"xmax"``, ``"zmin"``, ``"zmax"``.
                            At least one is required to ensure uniqueness of the solution.
    :type dirichlet_faces: list[str]
    :param diffusivity: Diffusivity values in m\\ :sup:`2`.s\\ :sup:`-1`. 
                        A list of two values is expected, diffusivity[0] is the value below the base level,
                        diffusivity[1] is the value above the base level. 
                        Default is :math:`10^{-10}` m\\ :sup:`2`.s\\ :sup:`-1` and :math:`10^{-6}` m\\ :sup:`2`.s\\ :sup:`-1`.
    :type diffusivity: list[float]
    :param baselevel: Base level in m. Default is 0.0.
    :param model_name: Model name.
    :type model_name: str
    :raises ValueError: If an unrecognized face is passed.
  """
  def __init__(self, dirichlet_faces:list[str],
               diffusivity:list[float]=[1.0e-10,1.0e-6],
               baselevel:float=0.0, 
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
    self.baselevel   = baselevel
    return
  
  def sprint_option(self):
    s = "########### Surface processes ###########\n"
    s += f"-{self.model_name}_{self.prefix}_apply_surface_diffusion # activate diffusion at surface\n"
    s += f"-{self.model_name}_{self.prefix}_varying_diffusivity\n"
    s += f"-{self.model_name}_{self.prefix}_diffusivity " + ','.join(map(str, self.diffusivity)) + '\n'
    s += f"-{self.model_name}_{self.prefix}_baselevel {self.baselevel:g}\n"
    s += "# spm diffusion boundary conditions minimum 1 must be passed\n"
    for face in self.bcs:
      s += f"-{self.model_name}_{self.prefix}_diffusion_dirichlet_{face}\n"
    return s
  
  def __str__(self):
    s = "Surface processes: Diffusion with base level\n"
    s += f"\tDiffusivity: {self.diffusivity}\n"
    s += f"\tBaselevel: {self.baselevel}\n"
    s += f"\tDirichlet faces: {self.bcs}\n"
    return s
  
def test():
  c = SPMDiffusionBaselvl(['xmin','xmax'],[1.0e-12,1.0e-5],10.0)
  s = c.sprint_option()
  print(s)

if __name__ == "__main__":
  test()