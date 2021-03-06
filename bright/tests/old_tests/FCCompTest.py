from __future__ import print_function
import sys
import os
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
    print("Do Calc Empty: " + str( fc.calc() ))
    testCD = {922350: 10.0, 10010: 1.0}
    print("Do Calc Dictonary: " + str( fc.calc(testCD) ))
    testMS = bright.MassStream({942390: 10.0, 80160: 20.0})
    print("Do Calc MassStream: " + str( fc.calc(testMS) ))
    return

fc0 = bright.FCComp()
print("Empty Fuel Cycle Component")
printFCComp(fc0)
print("")

print(bright.track_isos())
bright.track_isos( [922350, 942390, 10010] )
print(bright.track_isos())



fc1 = bright.FCComp()
print("Isotope track with no name...")
printFCComp(fc1)
print("")

bright.verbosity(1)
print(bright.verbosity())

fc2 = bright.FCComp("Isotope track with name!")
printFCComp(fc2)
print("")

bright.verbosity(100)
print(bright.verbosity())

ptrack = ["My first Param", "Another Param"]
fc3 = bright.FCComp(ptrack)
print("Full track with no name...")
printFCComp(fc3)
print("")

fc4 = bright.FCComp(ptrack, "Full track with name!")
printFCComp(fc4)
print("")


if os.name == "posix":
        subprocess.call("rm -r *Isos.txt *Params.txt", shell=True)
elif os.name == "nt":
        for f in os.listdir('.'):
                if "Isos.txt" in f:
                        os.remove(f)
                elif "Params.txt" in f:
                        os.remove(f)
