import numpy as np
from amuse.units import units
from amuse.units import quantities
from amuse.units import constants
from amuse.units import nbody_system
from amuse.ext.bridge import bridge
from amuse.community.ph4.interface import ph4
from amuse.community.bhtree.interface import BHTree


code_tree = BHTree()
code_direct = ph4()

print(code_tree)
print(code_direct)
