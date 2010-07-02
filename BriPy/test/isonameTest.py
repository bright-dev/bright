import os

import sys
if os.name == "posix":
    sys.path.append("../build/lib.linux-x86_64-2.6")
if os.name == "nt":
    sys.path.append("../build/lib.win32-2.6")

from isoname import *

print( LLAAAM_2_zzaaam("U234") )

print( mixed_2_zzaaam_List(["U235", 94239, 10010]) )
print( mixed_2_LLAAAM_List(["U235", 94239, 10010]) )
print( mixed_2_MCNP_List(["U235", 94239, 10010]) )
