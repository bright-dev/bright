from __future__ import print_function
import os
import sys
if os.name == "posix":
        sys.path.append("../build/lib.linux-x86_64-2.6")
elif os.name == "nt":
        sys.path.append("../build/lib.win32-2.6")

import subprocess

from bright import *

def printFCComp(fc):
    print("Name: " + fc.name)
    print("Params2Track: " + str(fc.track_params))
    print("ms_feed: " + str(fc.ms_feed))
    print("ms_prod: " + str(fc.ms_prod))
    print("params_prior_calc: " + str(fc.params_prior_calc))
    print("params_after_calc: " + str(fc.params_after_calc))
    print("Pass Number: " + str(fc.pass_num))
    return

def printReactorParameters(rp):
    print(rp.batches)
    print(rp.flux)
    print(rp.fuel_form)
    print(rp.coolant_form)
    print(rp.fuel_density)
    print(rp.coolant_density)
    print(rp.pnl)
    print(rp.BUt)
    print(rp.use_disadvantage_factor)
    print_ms()
    print(rp.radius)
    print(rp.pitch)
    print(rp.open_slots)
    print(rp.total_slots)
    print_ms()
    return

def printReactorVals(r):
    print("Batches", r.B)
    print("phi",r.phi)
    print("fuel_chemical_form", r.fuel_chemical_form)
    print("coolant_chemical_form", r.coolant_chemical_form)
    print("rhoF", r.rhoF)
    print("rhoC", r.rhoC)
    print("P_NL", r.P_NL)
    print("target_BU", r.target_BU)
    print("use_zeta",  r.use_zeta)
    print_ms()
    print("r", r.r)
    print("l", r.l)
    print("S_O", r.S_O)
    print("S_T", r.S_T)
    print("VF",  r.VF)
    print("VC",  r.VC)
    print_ms()
    return

track_isos([922350, 942390, 10010, 80160])

print("Fast Reactor Default Parameters:")
printReactorParameters(fr_defaults())

r0 = FastReactor1G()
print("Empty Reprocessing Cycle Component")
printFCComp(r0)
print("")

printReactorVals(r0)

r1 = FastReactor1G("FR.h5")
print("SE with no name...")
printFCComp(r1)
print("")

r2 = FastReactor1G("FR.h5", "Reactor2")
printFCComp(r2)
print("")


print(r2.F)
print_ms()
r2.loadlib("LWR.h5")
print(r2.F)
print_ms()

NewParams = fr_defaults()
NewParams.batches = 6
NewParams.fuel_density = 1200.0
NewParams.use_disadvantage_factor = True
NewParams.total_slots = 400

printReactorParameters(NewParams)

r3 = FastReactor1G(NewParams, "Reactor3")
printFCComp(r3)
printReactorVals(r3)
print_ms()

r4 = FastReactor1G("LWR.h5", NewParams, "Reactor4")
printFCComp(r4)
printReactorVals(r4)
print_ms()

if os.name == "posix":
        subprocess.call("rm -r *Isos.txt *Params.txt", shell=True)
elif os.name == "nt":
        for f in os.listdir('.'):
                if "Isos.txt" in f:
                        os.remove(f)
                elif "Params.txt" in f:
                        os.remove(f)
