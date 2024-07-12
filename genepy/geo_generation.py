import sys
from genepy import Domain
import numpy as np
import gmsh

class GeoPoints(Domain):
  def __init__(self,domain:Domain,mesh_size:float,name:str="gmsh_model") -> None:
    Domain.__init__(self,domain.dim,domain.O_num,domain.L_num,domain.n)
    self.mesh_size = mesh_size
    self.name      = name
    return

class GmshLayerBox(GeoPoints):
  def __init__(self, domain:Domain, mesh_size:float, interfaces, volumes, boundaries, name:str="gmsh_model") -> None:
    GeoPoints.__init__(self,domain,mesh_size,name)
    self.interfaces = np.asarray(interfaces, dtype=np.float64)
    self.volumes    = np.asarray(volumes,    dtype=np.int32)
    self.boundaries = np.asarray(boundaries, dtype=np.int32)
    self.n          = self.interfaces.shape[0]
    self.nlayer     = self.n-1
    self.nbcs       = self.boundaries.shape[0]

    if self.volumes.shape[0] != self.nlayer:
      raise ValueError(f"volumes must have one less element than interfaces, found {self.volumes.shape[0]} and {self.nlayer} respectively.")
    
    self.tags = {}
    self.tags["points"] = {"Oz": np.zeros(self.n, dtype=np.int32),
                           "Lz": np.zeros(self.n, dtype=np.int32)}
    self.tags["lines"] = {"horizontal": np.zeros(self.n, dtype=np.int32),
                          "vertical_Oz": np.zeros(self.nlayer, dtype=np.int32),
                          "vertical_Lz": np.zeros(self.nlayer, dtype=np.int32)}
    self.tags["surfaces"] = {"Ox": np.zeros(self.nlayer, dtype=np.int32),
                             "extruded": []}

    self.create_gmsh_mesh()
    return

  def create_horizontal_lines(self):
    """
    Create all the points of a cross section of normal x at Oz and Lz as well as 
    all the horizontal line between the points store their tags in a numpy array 
    """
    for i in range(self.n):
      self.tags["points"]["Oz"][i] = gmsh.model.geo.addPoint(self.O_num[0], self.interfaces[i], self.O_num[2], self.mesh_size)
      self.tags["points"]["Lz"][i] = gmsh.model.geo.addPoint(self.O_num[0], self.interfaces[i], self.L_num[2], self.mesh_size)
      self.tags["lines"]["horizontal"][i] = gmsh.model.geo.addLine(self.tags["points"]["Oz"][i], self.tags["points"]["Lz"][i])
    return

  def create_vertical_lines(self):
    """
    Create all the vertical lines between the points at Oz and Lz and store their tags.
    """
    for i in range(self.nlayer):
      self.tags["lines"]["vertical_Oz"][i]= gmsh.model.geo.addLine(self.tags["points"]["Oz"][i],
                                                                   self.tags["points"]["Oz"][i+1])
      self.tags["lines"]["vertical_Lz"][i]= gmsh.model.geo.addLine(self.tags["points"]["Lz"][i],
                                                                   self.tags["points"]["Lz"][i+1])
    return
  
  def create_surfaces(self):
    """
    Create all the close surface (rectangle) by joining the vertical and horizontal lines.
    """
    for i in range(self.nlayer):
      tag = gmsh.model.geo.addCurveLoop([self.tags["lines"]["horizontal"][i],
                                         self.tags["lines"]["vertical_Lz"][i],
                                         -self.tags["lines"]["horizontal"][i+1],
                                         -self.tags["lines"]["vertical_Oz"][i]])
      self.tags["surfaces"]["Ox"][i] = gmsh.model.geo.addPlaneSurface([tag])
    return
  
  def extrude_section(self):
    """
    Create all the volumes based on the cross section using gmsh extrude function
    store the extruded volume in the list ext_surf
    extruded volume is a list of list that contains the dimension [i][0] and the tag [i][1] for:

    .. code-block:: python

      [0][1] # --> tag of the surface of normal x+  that has been extruded
      [1][1] # --> tag of the volume 
      [2][1] # --> tag of the surface of normal y- 
      [3][1] # --> tag of the surface of normal z-
      [4][1] # --> tag of the surface of normal y+
      [5][1] # --> tag of the surface of normal z+
    
    """
    for i in range(self.nlayer):
      self.tags["surfaces"]["extruded"].append(gmsh.model.geo.extrude([(2,self.tags["surfaces"]["Ox"][i])],self.L_num[0]-self.O_num[0],0,0))
    gmsh.model.geo.synchronize()
    return

  def create_physical_volumes(self):
    """
    Using the user provided list of volumes, gather all the volumes of the same `pTatin3d`_ region in a 
    list and, for each region, add a physical group of dimension 3 gathering the list content.
    """
    dim = 3
    unique_values = np.unique(self.volumes)
    for tag in unique_values:
      tag_list = []
      for i in range(self.nlayer):
        if self.volumes[i] == tag:
          tag_list.append(self.tags["surfaces"]["extruded"][i][1][1])
      gmsh.model.addPhysicalGroup(dim, tag_list)
    gmsh.model.geo.synchronize()
    return
  
  def create_physical_surfaces(self):
    """ 
    Using the user provided boundaries list, gather all the faces with same `pTatin3d`_ bctype in a list 
    and add a physical group of dimension 2 for each list.
    Note the choice to separate each individual gmsh volume for the face of normal x and the face of normal z.
    """
    # Bottom boundary
    dim = 2
    gmsh.model.addPhysicalGroup(dim, [self.tags["surfaces"]["extruded"][0][2][1]])
    # Other boundaries
    unique_values = np.unique(self.boundaries)
    for tag in unique_values:
      tag_nx  = []
      tag_nz  = []
      for i in range(self.nbcs):
        if self.boundaries[i] == tag:
          tag_nx.append(self.tags["surfaces"]["Ox"][i])
          tag_nx.append(self.tags["surfaces"]["extruded"][i][0][1])
          tag_nz.append(self.tags["surfaces"]["extruded"][i][3][1])
          tag_nz.append(self.tags["surfaces"]["extruded"][i][5][1])
      gmsh.model.addPhysicalGroup(dim,tag_nx)
      gmsh.model.addPhysicalGroup(dim,tag_nz)
    gmsh.model.geo.synchronize()
    return
  
  def write_mesh(self):
    gmsh.model.mesh.generate()
    # generate a msh and a geo_unrolled file
    gmsh.write(self.name+".geo_unrolled")
    gmsh.write(self.name+".msh")
    # show the result in gmsh graphical interface
    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
    gmsh.finalize()
    return
  
  def create_gmsh_mesh(self):
    gmsh.initialize(sys.argv)
    gmsh.model.add(self.name)
    self.create_horizontal_lines()
    self.create_vertical_lines()
    self.create_surfaces()
    gmsh.model.geo.synchronize()
    self.extrude_section()
    self.create_physical_volumes()
    self.create_physical_surfaces()
    self.write_mesh()
    return

def test_box():
  # mesh size
  lc = 5e5
  # Domain
  Ox = 0. 
  Lx = 1500e3
  Oy = -650e3
  Ly = 0e3
  Oz = 0e3
  Lz = 1500e3
  domain = Domain(3, 
                  np.array([Ox,Oy,Oz],dtype=np.float64), 
                  np.array([Lx,Ly,Lz],dtype=np.float64), 
                  np.array([2,2,2],dtype=np.int32)
                 )
  # depth of interfaces
  interfaces = [Oy,-400e3,-120e3,-100e3,-80e3,-60e3,-40e3,-20e3,Ly]
  # tags of the volumes ordered by interfaces
  volumes    = [0,1,2,3,2,3,4,5]
  # tags of the boundaries
  bcs        = [0,0,1,1,1,1,2,2]

  box = GmshLayerBox(domain, lc, interfaces, volumes, bcs, name="multilyermsh")
  return

if __name__ == "__main__":
  test_box()