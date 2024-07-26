#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: mesh_refinement.py
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

import numpy as np
from genepy.initial_conditions import domain

class MeshRefinement(domain.Domain):
  """
  .. py:class:: MeshRefinement(Domain, refinement_params)

    Class to refine a mesh in one or more directions. 
    The refinement is done by linear interpolation of the coordinates of the mesh in the specified directions.
    The class is a subclass of :py:class:`Domain <genepy.initial_conditions.domain.Domain>`

    :param Domain Domain: domain to refine
    :param dict refinement_params: dictionary containing the refinement parameters. The dictionary must have the following structure:
    
    .. code-block:: python

      refinement_params = {"x": {"x_initial": np.array([...],dtype=np.float64), "x_refined": np.array([...],dtype=np.float64)},
                           "y": {"x_initial": np.array([...],dtype=np.float64), "x_refined": np.array([...],dtype=np.float64)},
                           "z": {"x_initial": np.array([...],dtype=np.float64), "x_refined": np.array([...],dtype=np.float64)}

    or

    .. code-block:: python

      refinement_params = {"x": {"mesh_fraction": np.array([...],dtype=np.float64), "x_refined": np.array([...],dtype=np.float64)}}

    The keys of the dictionary are the directions to refine (``"x"``, ``"y"``, ``"z"``) 
    and the values are dictionaries with the keys:
    
    1. ``"x_initial"`` and ``"x_refined"`` 
    2. ``"mesh_fraction"`` and ``"x_refined"``

    In the first case, the values of these keys are numpy arrays of the initial and refined coordinates in the specified direction. 
    In the second case, the values of these keys are numpy arrays of the mesh fractions and the refined coordinates in the specified direction.
    The arrays must have the same shape.

    Attributes
    ----------
    
    .. py:attribute:: params
      :type: dict
      :canonical: genepy.initial_conditions.mesh_refinement.MeshRefinement.params

      Refinement parameters

    .. py:attribute:: dirmap
      :type: dict
      :canonical: genepy.initial_conditions.mesh_refinement.MeshRefinement.dirmap

      Dictionary to map the direction to the index of the direction in the mesh i.e., ``"x":0``, ``"y":1``, ``"z":2``

    Methods
    -------
  """
  def __init__(self, Domain, refinement_params) -> None:
    domain.Domain.__init__(self,Domain.dim,Domain.O_num,Domain.L_num,Domain.n)
    self.dirmap = {"x":0,"y":1,"z":2}
    self.params = refinement_params
    for d in self.params:
      if d not in ["x","y","z"]:
        raise ValueError(f"Invalid direction, possible values are \"x\", \"y\", \"z\", found: {d}")
      for p in self.params[d]:
        if p not in ["x_initial","x_refined","mesh_fraction"]:
          raise ValueError(f"Invalid parameter, possible values are \"x_initial\", \"x_refined\", \"mesh_fraction\" found: {p}")

      if "x_initial" in self.params[d]:
        self.params[d]["x_initial"] = np.asarray(self.params[d]["x_initial"],dtype=np.float64)
      elif "mesh_fraction" in self.params[d]:
        self.params[d]["x_initial"] = np.cumsum(self.params[d]["mesh_fraction"]) * (self.L[self.dirmap[d]] - self.O[self.dirmap[d]]) + self.O[self.dirmap[d]]
      else:
        raise ValueError(f"Direction {d}: either \"x_initial\" or \"mesh_fraction\" must be provided")
      self.params[d]["x_refined"] = np.asarray(self.params[d]["x_refined"],dtype=np.float64)
      if self.params[d]["x_initial"].shape != self.params[d]["x_refined"].shape:
        raise ValueError(f"Direction {d}: array in \"x_initial\" of shape {self.params[d]['x_initial'].shape} and \"x_refined\" {self.params[d]['x_refined'].shape} must have the same shape")
    return

  def __str__(self) -> str:
    s  = f"{self.__class__.__name__}:\n"
    for d in self.params:
      s += f"\tDirection: {d}\n"
      for p in self.params[d]:
        s += f"\t\t{p}:\t{self.params[d][p]}\n"
        if p == "x_initial" or p == "x_refined":
          s += f"\t\t{p} normalized:\t{self.normalize(self.params[d][p],self.dirmap[d])}\n"
    return s
  
  def normalize(self,x,dim):
    """
    normalize(x,dim)
    Normalize the coordinates of the mesh in the specified direction

    :param np.ndarray x: coordinates to normalize
    :param int dim:      direction of the coordinates to normalize

    :return: normalized coordinates
    :rtype: np.ndarray
    """
    return (x - self.O[dim]) / (self.L[dim] - self.O[dim])
  
  def refine_direction(self,dim,x_initial,x_refined):
    """
    refine_direction(dim,x_initial,x_refined)
    Refine the mesh in the specified direction using linear interpolation.
    :attr:`num_coor <genepy.initial_conditions.domain.Domain.num_coor>` is updated with the refined coordinates.

    :param int dim:              direction of the mesh to refine (``0:x``, ``1:y``, ``2:z``)
    :param np.ndarray x_initial: initial coordinates 
    :param np.ndarray x_refined: refined coordinates
    """
    # convert tuple to list (tuples are immutable)
    coor = list(self.num_coor)
    # get all points in the dim direction
    xc_ini = coor[dim]
    # normalize
    xc_ini = self.normalize(xc_ini,dim)
    # linear interpolation
    xc_refined = np.interp(xc_ini,x_initial,x_refined)
    # rescale
    coor[dim] = self.O[dim] + xc_refined * (self.L[dim] - self.O[dim])
    # restore as a tuple
    self.num_coor = tuple(coor)
    return
  
  def refine(self):
    """
    refine(self)
    Refine the mesh in the specified directions in :attr:`params <MeshRefinement.params>`.
    """
    for d in self.params:
      x_initial = self.normalize(self.params[d]["x_initial"],self.dirmap[d])
      x_refined = self.normalize(self.params[d]["x_refined"],self.dirmap[d])
      self.refine_direction(self.dirmap[d],x_initial,x_refined)
    return
  
  def sprint_option(self,model_name:str):
    """
    sprint_option(self,model_name:str)
    Returns a string formatted for pTatin3d input file using `PETSc <https://petsc.org>`_ options format.
    
    :param str model_name: name of the model to include in the options

    :return: string with the options
    """
    prefix = "refinement"
    s  = f"###### Mesh {prefix} ######\n"
    s += f"-{model_name}_{prefix}_apply # activate mesh {prefix}\n"

    # count number of directions
    n = 0
    for d in self.params:
      n += 1
    s += f"-{model_name}_{prefix}_ndir {n} # number of directions (x,y,z) being refined\n"

    # write directions
    directions = list(self.params.keys())
    s += f"-{model_name}_{prefix}_dir "
    for i in range(n-1):
      s += f"{self.dirmap[directions[i]]},"
    s += f"{self.dirmap[directions[n-1]]} # directions of refinement (0:x, 1:y, 2:z)\n"

    # for each direction, report the number of points and the normalized xp,f(xp)
    for d in self.params:
      x_initial = self.normalize(self.params[d]["x_initial"],self.dirmap[d])
      x_refined = self.normalize(self.params[d]["x_refined"],self.dirmap[d])

      npoints = x_initial.shape[0]
      s += f"-{model_name}_{prefix}_npoints_{self.dirmap[d]} {npoints} # number of points for interpolation\n"
      s += f"-{model_name}_{prefix}_xref_{self.dirmap[d]} "
      for i in range(npoints-1):
        s += f"{x_initial[i]},"
      s += f"{x_initial[npoints-1]} # xp\n"
      s += f"-{model_name}_{prefix}_xnat_{self.dirmap[d]} "
      for i in range(npoints-1):
        s += f"{x_refined[i]},"
      s += f"{x_refined[npoints-1]} # f(xp)\n"
    return s

def test():
  from genepy import writers
  O = np.array([0,-250,0],dtype=np.float64)
  L = np.array([600,0,600],dtype=np.float64)
  n = np.array([45,18,33],dtype=np.int32)
  d = domain.Domain(3,O,L,n)
  #w = writers.WriteVTS(d,vtk_fname="mesh_refinement.vts")

  #"x": {"x_initial": np.array([0,12,60,540,588,600], dtype=np.float64),
  #      "x_refined": np.array([0,100,150,450,500,600], dtype=np.float64)},
  refinement = {"x": {"mesh_fraction": np.array([0,0.02,0.08,0.8,0.08,0.02], dtype=np.float64),
                      "x_refined":     np.array([0,100,150,450,500,600], dtype=np.float64)},
                "y": {"x_initial":     np.array([-250,-180,-87.5,0], dtype=np.float64),
                      "x_refined":     np.array([-250,-50,-16.25,0], dtype=np.float64)},
                "z": {"x_initial":     np.array([0,12,60,540,588,600], dtype=np.float64),
                      "x_refined":     np.array([0,180,240,360,420,600], dtype=np.float64)}}

  m = MeshRefinement(d,refinement)
  print(m)
  options = m.sprint_option("model_GENE3D")
  print(options)
  m.refine()
  w = writers.WriteVTS(m,vtk_fname="mesh_refinement.vts")
  w.write_vts()

if __name__ == "__main__":
  test()
