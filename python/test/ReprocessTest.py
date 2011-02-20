from __future__ import print_function
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
    print("Do Calc Empty: " + str( fc.doCalc() ))
    testCD = {922350: 10.0, 10010: 1.0}
    print("Do Calc Dictonary: " + str( fc.doCalc(testCD) ))
    testMS = bright.MassStream({942390: 10.0, 80160: 20.0})
    print("Do Calc MassStream: " + str( fc.doCalc(testMS) ))
    return

rp0 = bright.Reprocess()
print("Empty Reprocessing Cycle Component")
printFCComp(rp0)
print("")

SE1 = {92: 0.99, 94: 0.9}
bright.isos2track([922350, 942390, 10010])

rp1 = bright.Reprocess({})
rp1.initialize(SE1)
print("SE with no name...")
printFCComp(rp1)
print("")

SE2 = {922350: 0.99, 942390: 0.9}
rp2 = bright.Reprocess({}, "SE with name!")
rp2.initialize(SE2)
printFCComp(rp2)
print("")

SE3 = {"U": 0.99, "PU239": 0.9}
rp3 = bright.Reprocess(SE3)
print("String SE with no name...")
printFCComp(rp3)
print("")

SE4 = {"U": 0.99, "PU239": 0.9, "80160": 1.0}
rp4 = bright.Reprocess(SE3, "String SE with name!")
printFCComp(rp4)
print("")

print("Starting Real Reprocessing trial...")
ms = bright.MassStream("MassStreamtry02.txt", -1, "RealMS")

isotrack = ms.comp.keys()
SE = {"U": 0.99, "PU": 0.9, "CM": 0.8}
R = bright.Reprocess(SE, "RealRep")
print(str(R.sepeff))
print("")
Rsef = R.sepeff
Rsef[952421] = 0.5
R.sepeff = Rsef
print(str(R.sepeff))
print("")

R.doCalc(ms)
print(str(R.IsosIn))
print(str(R.IsosOut))
print("")

print(str(R.ParamsIn))
print(str(R.ParamsOut))
print("")
R.setParams()
print(str(R.ParamsIn))
print(str(R.ParamsOut))
print("")


if os.name == "posix":
        subprocess.call("rm -r *Isos.txt *Params.txt", shell=True)
elif os.name == "nt":
        for f in os.listdir('.'):
                if "Isos.txt" in f:
                        os.remove(f)
                elif "Params.txt" in f:
                        os.remove(f)
