import os

version = [1,1,3]

from .domain          import *
from .rotation        import *
from .bcs             import *
from .utils           import *
from .gaussian        import *
from .writers         import *
from .bc_inversion    import *
from .mesh_refinement import *

# Material parameters
from .material_params.materials  import *
from .material_params.density    import *
from .material_params.plasticity import *
from .material_params.softening  import *
from .material_params.viscosity  import *
from .material_params.energy     import *
from .material_params.regions    import *