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
    print("Params2Track: " + str(fc.track_params))
    print("ms_feed: " + str(fc.ms_feed))
    print("ms_prod: " + str(fc.ms_prod))
    print("params_prior_calc: " + str(fc.params_prior_calc))
    print("params_after_calc: " + str(fc.params_after_calc))
    print("Pass Number: " + str(fc.pass_num))
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

st1.ms_feed = ms
print("Empty Decay")
st1.calc()
print(str( st1.ms_prod ))

st1.calc(cd)
print("CompDict Decay")
print(str( st1.ms_prod ) )

st1.decay_time =  st1.decay_time * 10000
st1.calc(ms)
print("MassStream Decay")
print(str( st1.ms_prod ) )

st2.ms_feed = ms 
t = 3600 * 24 * 10000000
st2.calc(t)
print("Double Decay")
print(str( st2.ms_prod ) )

st2.calc(cd, t)
print("CompDict, Double Decay")
print(str( st2.ms_prod ) )

st2.calc(ms, t*10000)
print("MassStream, Double Decay")
print(str( st2.ms_prod ) )


if os.name == "posix":
        subprocess.call("rm -r *Isos.txt *Params.txt", shell=True)
elif os.name == "nt":
        for f in os.listdir('.'):
                if "Isos.txt" in f:
                        os.remove(f)
                elif "Params.txt" in f:
                        os.remove(f)
