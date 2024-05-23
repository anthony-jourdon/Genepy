import os

version = [1,1,3]

from .initial_conditions.domain          import *
from .initial_conditions.gaussian        import *
from .initial_conditions.mesh_refinement import *
from .initial_conditions.ics             import *

from .rotation        import *
from .utils           import *
from .writers         import *
from .bc_inversion    import *

# Boundary conditions
from .boundary_conditions.velocity   import *
from .boundary_conditions.bcs        import *
from .boundary_conditions.dirichlet  import *
from .boundary_conditions.navierslip import *

# Markers
from .markers.popctrl import *
from .markers.pswarm  import *

# Material parameters
from .material_params.materials  import *
from .material_params.density    import *
from .material_params.plasticity import *
from .material_params.softening  import *
from .material_params.viscosity  import *
from .material_params.energy     import *
from .material_params.regions    import *

# Options
from .ptatin_options    import *
from .surface_processes import *