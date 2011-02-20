import os
import sys
if os.name == "posix":
        sys.path.append("../build/lib.linux-x86_64-2.6")
elif os.name == "nt":
        sys.path.append("../build/lib.win32-2.6")

import subprocess
import bright

def printFCComp(fc):
    print("Name: " + fc.name)
    print("Params2Track: " + str(fc.params2track))
    print("IsosIn: " + str(fc.IsosIn))
    print("IsosOut: " + str(fc.IsosOut))
    print("ParamsIn: " + str(fc.ParamsIn))
    print("ParamsOut: " + str(fc.ParamsOut))
    print("Pass Number: " + str(fc.PassNum))
    return

st0 = bright.Storage()
print("Empty Storage Component")
printFCComp(st0)
print("")

bright.track_isos([922350, 942390, 10010])
cd = {922350: 10.0, 10010: 1.0}
ms = bright.MassStream("MassStreamtry02.txt")

st1 = bright.Storage()
print("Storage No Name...")
printFCComp(st1)
print("")

bright.track_isos(ms.comp.keys())
st2 = bright.Storage("Storage With Name!")
printFCComp(st2)
print("")

print("ST1 Decay Time: " + str(st1.decay_time))
st1.decay_time = 3600 * 24 * 10
print("ST1 Decay Time: " + str(st1.decay_time))
print("")

st1.IsosIn = ms
print("Empty Decay")
st1.doCalc()
print(str( st1.IsosOut ))

st1.doCalc(cd)
print("CompDict Decay")
print(str( st1.IsosOut ) )

st1.decay_time =  st1.decay_time * 10000
st1.doCalc(ms)
print("MassStream Decay")
print(str( st1.IsosOut ) )

st2.IsosIn = ms 
t = 3600 * 24 * 10000000
st2.doCalc(t)
print("Double Decay")
print(str( st2.IsosOut ) )

st2.doCalc(cd, t)
print("CompDict, Double Decay")
print(str( st2.IsosOut ) )

st2.doCalc(ms, t*10000)
print("MassStream, Double Decay")
print(str( st2.IsosOut ) )


if os.name == "posix":
        subprocess.call("rm -r *Isos.txt *Params.txt", shell=True)
elif os.name == "nt":
        for f in os.listdir('.'):
                if "Isos.txt" in f:
                        os.remove(f)
                elif "Params.txt" in f:
                        os.remove(f)
