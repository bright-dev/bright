#! /usr/bin/env python

import os
import sys
if os.name == "posix":
        sys.path.append("../../build/lib.linux-x86_64-2.6")
elif os.name == "nt":
        sys.path.append("../../build/lib.win32-2.6")
        os.putenv("HDF5_DISABLE_VERSION_CHECK", "2")

import shutil
from optparse import OptionParser

def main():
    parser = OptionParser()
    parser.add_option("-c", "--case", dest="case", help="Benchmark case to run.")
    parser.add_option("-p", "--calibrate", action="store_true", dest="calibrate", default=False, help="Calibrate non-leakage probability.")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="prints extra information.")
    options, args = parser.parse_args()

    if options.case == "LWR_NEA":
        name = options.case
    elif options.case == "LWR_VIS51":
        name = options.case
    else:
        print "Case not valid, please pick from: LWR_NEA, or LWR_VIS51."
        raise SystemExit

    shutil.copy(name + "_fcparams.py", "fcparams.py")

    import bright
    import tables

    from fcparams import lwr_params
    from fcparams import Quiet

    #Various Variables
    snf_need = []
    if (not Quiet) or (options.verbose):
        bright.verbosity(100)


    #redefine isotrak
    trackfile = tables.openFile("../LWR.h5", 'r')
    itrack = trackfile.root.ToIso_zz.read()
    trackfile.close()
    bright.isos2track(itrack)

    if name == "LWR_NEA":
        #NEA
        U234 = bright.MassStream({922340: 1.0}, 0.00032, "U234")
        U235 = bright.MassStream({922350: 1.0}, 0.03600, "U235")
        U236 = bright.MassStream({922350: 1.0}, 0.00016, "U235")
        U238 = bright.MassStream({922380: 1.0}, 0.96352, "U238")
    elif name == "LWR_VIS51":
        #VISION
        U234 = bright.MassStream({922340: 1.0}, 3.439849E-04, "U234")
        U235 = bright.MassStream({922350: 1.0}, 4.299811E-02, "U235")
        U236 = bright.MassStream({922350: 1.0}, 0.000000E+00, "U235")
        U238 = bright.MassStream({922380: 1.0}, 9.566579E-01, "U238")
    else:
        print "Case not valid, please pick from: LWR_NEA, or LWR_VIS51."
        raise SystemExit

    #######################
    ### LWR Computation ###
    #######################

    #Fuel Cycle Components
    LWR = bright.LightWaterReactor1G("../LWR.h5", lwr_params, name)

    def LWR_delR_BU_(ms):
        "Calculates the delta Reaction Rates at the target burnup."
        LWR.IsosIn = ms
        LWR.foldMassWeights()
        dR = LWR.batchAve(lwr_params.BUt, "p") - LWR.batchAve(lwr_params.BUt, "d")
        return dR

    def Run_PNL(temp_pnl):
        LWR.P_NL = temp_pnl

        delR_U235 = LWR_delR_BU_(U235)
        delR_U238 = LWR_delR_BU_(U238)

        #Calculate delta R for the Guess
        LWR_CoreInput = U238 + U235 + U234 + U236
        LWR_CoreInput.name = "LWR_CoreInput"
        LWR_CoreInput.normalize()
        LWR_delR_Guess = LWR_delR_BU_(LWR_CoreInput)

        k = LWR.batchAveK(lwr_params.BUt)
        n = 0
        if not Quiet:
            print str(1) + ")",  k, 

        while 0.001 < abs(1.0 - k) and n < 10:
            #Adjust Masses based on pertubation guess.
            LWR_DeltaM_U238 = - LWR_delR_Guess / (delR_U238 - delR_U235)
            U238.mass = U238.mass + LWR_DeltaM_U238
            U235.mass = U235.mass - LWR_DeltaM_U238

            #Recalculate core parameters for new masses guess
            LWR_CoreInput = U238 + U235 + U234 + U236
            LWR_CoreInput.name = "LWR_CoreInput"
            LWR_delR_Guess = LWR_delR_BU_(LWR_CoreInput)
            k = LWR.batchAveK(lwr_params.BUt)
            n = n+1
            if not Quiet:
                print k, 
        if not Quiet:
            print
            print

        #Calculate and write output
        LWR.BUd_BisectionMethod()
        LWR.calcOutIso()
        LWR.writeout()

    if options.calibrate:
        LWR_CoreInput = U238 + U235 + U234 + U236
        LWR_CoreInput.name = "LWR_CoreInput"
        LWR_CoreInput.normalize()

        LWR.IsosIn = LWR_CoreInput

        print(LWR.IsosIn)

        LWR.Calibrate_PNL_2_BUd()

        print
        print "Non-Leakage Probability = ", LWR.P_NL

        print
        Run_PNL(LWR.P_NL)

        with open(name + "_pnl.txt", 'w') as f:
            f.write(str(LWR.P_NL))
    else:
        print
        Run_PNL(0.98)


if __name__ == "__main__":
    main()

"""print
print "SigmaFa = ", LWR.SigmaFa_F_
print
print "SigmaCa = ", LWR.SigmaCa_F_
print
print "SigmaFtr = ", LWR.SigmaFtr_F_
print
print "SigmaCtr = ", LWR.SigmaCtr_F_
print
print "kappaF = ", LWR.kappaF_F_
print
print "kappaC = ", LWR.kappaC_F_
print
print "zeta = ", LWR.zeta_F_
"""

