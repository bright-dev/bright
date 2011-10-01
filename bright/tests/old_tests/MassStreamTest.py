from __future__ import print_function
import os
import sys

if os.name == "posix":
    sys.path.append("../build/lib.linux-x86_64-2.6/")
elif os.name == "nt":
    sys.path.append("../build/lib.win32-2.6")

import MassStream 

ms = MassStream.MassStream("MassStreamtry01.txt", 42.0, "Trial Stream")
msstr = str( ms )
print(msstr)

ms2 = MassStream.MassStream("MassStreamtry02.txt", 42.0, "Trial Stream2")
ms2.print_ms()
print(ms2)
print(ms2.mult_by_mass())
ms2.normalize()
print(ms2)

mdict = {10010: 10.0, 80160: 5.0}
ms3 = MassStream.MassStream(mdict, 15.0, "Trial")
ms3.print_ms()

ms4 = ms3.get_sub_streamStr(["H"], "H")
ms4.print_ms()

ms5 = ms3.get_sub_streamInt([8, 10])
ms5.print_ms()
print("")

print("Get U?")
ms6 = ms2.get_u("Uranium From ms2")
print("Got U!")
ms6.print_ms()

print("Get FP?")
ms7 = ms3.get_fp()
print("Got FP!")
print(ms7)
print("")

print("Operator Overloading!")
print(ms3 + 30.0)
print(90 + ms3)
print(ms3 + ms6)

print(ms3 * 2)
print(150 * ms3)
print(ms3 / 10)
print("")
