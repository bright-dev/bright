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
    print("Params2Track: " + str(fc.params2track))
    print("IsosIn: " + str(fc.IsosIn))
    print("IsosOut: " + str(fc.IsosOut))
    print("ParamsIn: " + str(fc.ParamsIn))
    print("ParamsOut: " + str(fc.ParamsOut))
    print("Pass Number: " + str(fc.PassNum))
    return

def printReactorParameters(rp):
    print(rp.batches)
    print(rp.flux)
    print(rp.FuelForm)
    print(rp.CoolantForm)
    print(rp.FuelDensity)
    print(rp.CoolantDensity)
    print(rp.pnl)
    print(rp.BUt)
    print(rp.useDisadvantage)
    print_ms()
    print(rp.Radius)
    print(rp.Length)
    print(rp.open_slots)
    print(rp.total_slots)
    print_ms()
    return

def printReactorVals(r):
    print("Batches", r.B)
    print("phi",r.phi)
    print("FuelChemicalForm", r.FuelChemicalForm)
    print("CoolantChemicalForm", r.CoolantChemicalForm)
    print("rhoF", r.rhoF)
    print("rhoC", r.rhoC)
    print("P_NL", r.P_NL)
    print("TargetBU", r.TargetBU)
    print("useZeta",  r.useZeta)
    print_ms()
    print("r", r.r)
    print("l", r.l)
    print("S_O", r.S_O)
    print("S_T", r.S_T)
    print("VF",  r.VF)
    print("VC",  r.VC)
    print_ms()
    return

ms = MassStream("MassStreamtry02.txt")
isos2track([922350, 942390, 10010, 80160, 50100, 50110])

print("Light Water Reactor Default Parameters:")
printReactorParameters(LWRDefaults())

r0 = LightWaterReactor1G()
print("Empty Reprocessing Cycle Component")
printFCComp(r0)
print("")

printReactorVals(r0)

r1 = LightWaterReactor1G("LWR.h5")
print("SE with no name...")
printFCComp(r1)
print("")

r2 = LightWaterReactor1G("LWR.h5", "Reactor2")
printFCComp(r2)
print("")


print(r2.F)
print_ms()
r2.loadLib("FR.h5")
print(r2.F)
print_ms()

NewParams = LWRDefaults()
#NewParams.batches = 6
#NewParams.FuelDensity = 1200.0
#NewParams.useDisadvantage = True
#NewParams.total_slots = 400
#NewParams.CoolantForm    = {'O16': 1.0, 'H1': 2.0}
NewParams.CoolantDensity = 1.0

printReactorParameters(NewParams)

r3 = LightWaterReactor1G(NewParams, "Reactor3")
printFCComp(r3)
printReactorVals(r3)
print_ms()

ms = MassStream("MassStreamtry03.txt")
isos2track([ 922350, 922380, 942390, 10010, 80160, 50100, 50110])

r4 = LightWaterReactor1G("LWR.h5", NewParams, "Reactor4")
printFCComp(r4)
printReactorVals(r4)
print_ms()

r4.IsosIn = ms
r4.foldMassWeights()
print(r4.MWF)
print(r4.MWC)

print_ms()
print("SigmaFa = ", r4.SigmaFa_F_)
print_ms()
print("SigmaCa = ", r4.SigmaCa_F_)
print_ms()
print("kappaF = ", r4.kappaF_F_)
print_ms()
print("kappaC = ", r4.kappaC_F_)
print_ms()
print("LatticeF = ", r4.LatticeF_F_)
print_ms()
print("LatticeE = ", r4.LatticeE_F_)
print_ms()
print("zeta = ", r4.zeta_F_)

if os.name == "posix":
        subprocess.call("rm -r *Isos.txt *Params.txt", shell=True)
elif os.name == "nt":
        for f in os.listdir('.'):
                if "Isos.txt" in f:
                        os.remove(f)
                elif "Params.txt" in f:
                        os.remove(f)
