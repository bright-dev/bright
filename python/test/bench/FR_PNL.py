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
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="Prints extra information.")
    options, args = parser.parse_args()

    if options.case == "FR_NEA":
        name = options.case
        bud_t = 27.3969072997
    elif options.case == "FR_VISp1":
        name = options.case
        bud_t = 176.6
    elif options.case == "FR_VISp5":
        name = options.case
        bud_t = 176.6
    else:
        print "Case not valid, please pick from: FR_NEA, FR_VISp1, or FR_VISp5."
        raise SystemExit

    shutil.copy(name + "_fcparams.py", "fcparams.py")

    import BriPy
    import tables

    from fcparams import fr_params
    from fcparams import Quiet

    #Various Variables
    snf_need = []
    if (not Quiet) or (options.verbose):
        BriPy.verbosity(100)

    #redefine isotrak
    trackfile = tables.openFile("../FR.h5", 'r')
    itrack = trackfile.root.ToIso_zz.read()
    trackfile.close()
    BriPy.isos2track(itrack)

    ######################
    ### FR Computation ###
    ######################
    InStream = BriPy.MassStream(name + '_Benchmark_In.txt', 1.0, "InStream")


    #Fuel Cycle Components
    FR = BriPy.FastReactor1G("../FR.h5", fr_params, name)

    def Run_PNL(temp_pnl):
        FR.P_NL = temp_pnl

        #Calculate output
        FR.IsosIn = InStream
        FR.foldMassWeights()
        FR.BUd_BisectionMethod()


    #Calibration proceeds by bisection method...
    pnl_a = 0.6
    Run_PNL(pnl_a)
    bud_a = FR.BUd
    sign_a = (bud_a - bud_t) / abs(bud_a - bud_t)

    pnl_b = 0.7
    Run_PNL(pnl_b)
    bud_b = FR.BUd
    sign_b = (bud_b - bud_t) / abs(bud_b - bud_t)

    DoA = 10.0**(-15)        #Degree of accuracy to carry out calculations to.
    q = 0
    while (DoA < abs(pnl_a - pnl_b)) and (DoA < abs(bud_a - bud_b)) and q < 100:
        pnl_c = (pnl_a + pnl_b) / 2.0
        Run_PNL(pnl_c)
        bud_c = FR.BUd
        sign_c = (bud_c - bud_t) / abs(bud_c - bud_t)

        q = q + 1

        if (sign_a == sign_c) and not (sign_b == sign_c):
            pnl_a = pnl_c 
            bud_a = bud_c
            sign_a = sign_c
        elif (sign_b == sign_c) and not (sign_a == sign_c):
            pnl_b = pnl_c 
            bud_b = bud_c
            sign_b = sign_c
        else:
            if not Quiet:
                print
                print "SOMEWHERE WHILE FINDING k SOMETHING WENT WRONG!!!"
                print "Here is some information that might help you debug ^_^"
                print "pnl_%(ltr)s = %(pnl)f\tBUd_%(ltr)s = %(bud)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'a', 'pnl': pnl_a, 'bud': bud_a, 'sign': sign_a}
                print "pnl_%(ltr)s = %(pnl)f\tBUd_%(ltr)s = %(bud)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'b', 'pnl': pnl_b, 'bud': bud_b, 'sign': sign_b}
                print "pnl_%(ltr)s = %(pnl)f\tBUd_%(ltr)s = %(bud)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'c', 'pnl': pnl_c, 'bud': bud_c, 'sign': sign_c}
                print

    if not Quiet:
        print
        print "Final Result of Burnup Bisection Method Calculation:"
        print "q = ", q
        print "pnl_%(ltr)s = %(pnl).16f\tBUd_%(ltr)s = %(bud)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'a', 'pnl': pnl_a, 'bud': bud_a, 'sign': sign_a}
        print "pnl_%(ltr)s = %(pnl).16f\tBUd_%(ltr)s = %(bud)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'b', 'pnl': pnl_b, 'bud': bud_b, 'sign': sign_b}
        print "pnl_%(ltr)s = %(pnl).16f\tBUd_%(ltr)s = %(bud)f\tsign_%(ltr)s = %(sign)f"%{'ltr': 'c', 'pnl': pnl_c, 'bud': bud_c, 'sign': sign_c}

    #Write output
    FR.calcOutIso()
    FR.writeout()

    with open(name + "_pnl.txt", 'w') as f:
        f.write(str(FR.P_NL))

    with open(name + "_trucr.txt", 'w') as f:
        f.write(str(FR.TruCR))

if __name__ == "__main__":
    main()
