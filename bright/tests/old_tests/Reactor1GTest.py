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

r0 = Reactor1G()
print("Empty Reprocessing Cycle Component")
printFCComp(r0)
print("")

track_isos([922350, 942390, 10010, 80160])

r1 = Reactor1G()
print("SE with no name...")
printFCComp(r1)
print("")

r2 = Reactor1G("SE with name!")
printFCComp(r2)
print("")

r2.loadlib("FR.h5")
#print(r2.F[0])
#print(r2.BUi_F_[10010])
#print(r2.pi_F_[922350])
#print(r2.di_F_[922350])
print(r2.Tij_F_.keys())
print_ms()

rp = ReactorParameters()
rp.batches = 3
rp.flux = 5.0**15
rp.fuel_form = {"U235": 1.0, "O16": 2.0}
rp.coolant_form = {"H1": 2.0, "O16": 1.0}
rp.fuel_density = 10.7
rp.coolant_density = 0.73
rp.pnl = 0.98
rp.BUt = 1.4
rp.use_disadvantage_factor = False

rp.radius = 0.411
rp.pitch = 1.4
rp.open_slots = 32
rp.total_slots = 156

r2.initialize(rp)

print("Batches", r2.B)
print("phi",r2.phi)
print("fuel_chemical_form", r2.fuel_chemical_form)
print("coolant_chemical_form", r2.coolant_chemical_form)
print("rhoF", r2.rhoF)
print("rhoC", r2.rhoC)
print("P_NL", r2.P_NL)
print("target_BU", r2.target_BU)
print("use_zeta",  r2.use_zeta)
print_ms()
print("r", r2.r)
print("l", r2.l)
print("S_O", r2.S_O)
print("S_T", r2.S_T)
print("VF",  r2.VF)
print("VC",  r2.VC)
print_ms()

r2.fold_mass_weights()
#print(r2.P_F_)

#r2.calc_Mj_F_()
#print(r2.Mj_F_)
#r2.calc_Mj_Fd_()

#print(r2.ms_prod)
r2.calc_ms_prod()
#print(r2.ms_prod)

r2.calcSubStreams()
#print(r2.ms_feed_u)
#print(r2.ms_feed_tru)
#print(r2.ms_feed_lan)
#print(r2.ms_feed_act)
#print(r2.ms_prod_u)
#print(r2.ms_prod_tru)
#print(r2.ms_prod_lan)
#print(r2.ms_prod_act)

r2.calc_tru_cr()

fp = r2.fluence_at_BU(10.0)

verbosity(100)
bAK = r2.batch_average(40.0, "P")
#print(bAK)

#print(r2.P_NL)
#r2.run_P_NL(0.5)
#print(r2.P_NL)



testCD = {922350: 10.0, 10010: 1.0}
print("Do Calc Dictonary: " + str( r2.calc(testCD) ))
testMS = MassStream({942390: 10.0, 80160: 20.0})
print("Do Calc MassStream: " + str( r2.calc(testMS) ))
print("Do Calc Empty: " + str( r2.calc() ))


if os.name == "posix":
        subprocess.call("rm -r *Isos.txt *Params.txt", shell=True)
elif os.name == "nt":
        for f in os.listdir('.'):
                if "Isos.txt" in f:
                        os.remove(f)
                elif "Params.txt" in f:
                        os.remove(f)
