#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  bc-pre-processing
#  filename: ics.py
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

from bcpy import Domain

class InitialConditions:
  """
  .. class:: InitialConditions(Domain, velocity, model_name="model_GENE3D", **kwargs)

    Abstract class to generate options for 
    `pTatin3d`_ 
    model initial conditions.

    :param Domain Domain: instance of the Domain class
    :param velocity: velocity function (str or sympy expression)
    :param str model_name: name of the model (default: model_GENE3D)
    :param kwargs: keyword arguments

    Attributes
    ----------

    .. py:attribute:: model_name
      :type: str
      :canonical: bcpy.initial_conditions.ics.InitialConditions.model_name

        Name of the model to include in the options

    .. py:attribute:: Domain
      :type: Domain
      :canonical: bcpy.initial_conditions.ics.InitialConditions.Domain

        Instance of the Domain class

    .. py:attribute:: u
      :type: velocity
      :canonical: bcpy.initial_conditions.ics.InitialConditions.u

        Velocity function (str or sympy expression)

    .. py:attribute:: kwargs
      :type: dict
      :canonical: bcpy.initial_conditions.ics.InitialConditions.kwargs

        Keyword arguments

    .. py:attribute:: possible_kwargs
      :type: list
      :canonical: bcpy.initial_conditions.ics.InitialConditions.possible_kwargs

        Possible keyword arguments

        - **mesh_refinement**: instance of the MeshRefinement class
        - **initial_strain**: instance of the Gaussian class

    Methods:
    --------
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
    """
    sprint_option(self)
    Returns a string formatted for `pTatin3d`_ input file using `PETSc`_ options format
    for all the initial conditions.

    :return: string formatted for `pTatin3d`_ input file
    :rtype: str
    """
    s = self.Domain.sprint_option(self.model_name)
    s += self.sprint_initial_velocity()
    for kw in self.possible_kwargs:
      if kw in self.kwargs:
        s += self.kwargs[kw].sprint_option(self.model_name)
    return s
  
  def sprint_initial_velocity(self) -> str:
    """
    sprint_initial_velocity(self)
    Returns a string formatted for `pTatin3d`_ input file using `PETSc`_ options format
    for the initial velocity field.

    :return: string formatted for `pTatin3d`_ input file
    :rtype: str
    """
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
