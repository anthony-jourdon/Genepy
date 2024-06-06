#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: materials.py
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

class MaterialConstants:
  """
  .. py:class:: MaterialConstants(model_name:str="model_GENE3D", region:int=0)

    Abstract class to generate the options for a given type and region.
    All the material parameters classes inherit from this class and use its
    :meth:`MaterialConstants.sprint_option` method to generate the options for the given material parameter.
    Thus, the classes attributes must be named as the options of the material parameter.
    If an instance class attribute is not an option of the material parameter, it must be
    excluded from the :meth:`MaterialConstants.sprint_option` method.

    :param str model_name: Name of the model, default is "model_GENE3D"
    :param int region: Region number, default is 0.

    .. note:: 
      Except if the :class:`MaterialConstants` class (or one of its children) 
      is used outside of the context of the :class:`Regions` class, the name of the model
      and the region should not be defined when calling the class.

    Attributes
    ----------

    .. py:attribute:: model_name
      :type: str

        Name of the model, default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number, default is 0
    
    Methods
    -------
  """
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    self.model_name = model_name
    self.region     = region
  
  def sprint_option(self):
    """
    sprint_option(self)
    Return the options of the material parameter as a string.
    The option generation is automatised based on the attributes of the class (and its children).

    :return: The options of the material parameter as a string.
    :rtype: str
    """
    attributes = vars(self)
    s = ""
    for p in attributes:
      if p in ['model_name','region']:
        continue
      if type(attributes[p]) is str: fmt = f"{attributes[p]}"
      else:                          fmt = f"{attributes[p]:g}"
      s += f"-{self.model_name}_{p}_{self.region} {fmt}\n"
    return s
