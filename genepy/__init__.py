#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: __init__.py
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
from .boundary_conditions.neumann    import *

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
