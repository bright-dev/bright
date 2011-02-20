#! /usr/bin/env python

import math
import shutil
import subprocess
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-c", "--case", dest="case", help="Benchmark case to run.")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="Prints extra information.")
options, args = parser.parse_args()

if options.case == "LWR_NEA":
    name = options.case
    type = "LWR"
    header = "LWR NEA/OECD Benchmark\n    Target BU: 40.0\n    Batches: 4\n"
elif options.case == "LWR_VIS51":
    name = options.case
    type = "LWR"
    header = "LWR Vision Benchmark\n    Target BU: 51.0\n    Batches: 3\n"
elif options.case == "FR_NEA":
    name = options.case
    type = "FR"
    header = "FR NEA/OECD Benchmark, Used NEA Physics of Plutonium vol 5\n    Target BU: 27.3969072997\n    Batches: 1\n    TRU_CR: 0.5\n"
elif options.case == "FR_VISp1":
    name = options.case
    type = "FR"
    header = "FR VISION First Pass Benchmark\n    Target BU: 176.6\n    Batches: 3\n    TRU_CR: 0.5\n"
elif options.case == "FR_VISp5":
    name = options.case
    type = "FR"
    header = "FR VISION Fifth Pass Benchmark\n    Target BU: 176.6\n    Batches: 3\n    TRU_CR: 0.5\n"
else:
    print "Case not valid, please pick from: LWR_NEA, LWR_VIS51, FR_NEA, FR_VISp1, or FR_VISp5."
    raise SystemExit

def GrabN(n):
    "Grabs the index n from both files."
    "Returns dictionary of tuples."

    dat_dict = {}
    dat_ord = []

    bench_f = open(name + "_Benchmark.txt")
    for bench_line in bench_f:
        bls = bench_line.split()
        if bls[0] == "Isotope":
            continue
        model_f = open('%sIsos.txt'%name, 'r')
        for model_line in model_f:
            mls = model_line.split()
            if bls[0] == mls[0]:
                dat_dict[mls[0]] = (float(bls[n]), float(mls[n]))
                dat_ord.append(mls[0])
                break
        model_f.close()
    bench_f.close()
    return dat_dict, dat_ord

def SummaryTable(name, n):
    dat, ord = GrabN(n)
    s = ""
    s = s + name + "\n"
    s = s + "{0:10}{1:^15}{2:^15}{3:^15}{4:^15}{5:^15}{6:^15}\n".format('Isotope', 'Benchmark', 'Model', 'Delta', 'Mean', 'Sigma', 'Fractional Dev')
    for iso in ord:
        if (dat[iso][0] == 0.0) and (dat[iso][1] == 0.0):
            continue

        delta =  dat[iso][0] - dat[iso][1]
        mean =  (dat[iso][0] + dat[iso][1]) / 2.0
        sigma = math.sqrt(0.5 * ( (dat[iso][0] - mean)**2 + (dat[iso][1] - mean)**2 ) )
        if mean == 0.0:
            frac_dev = 0.0
        else:
            frac_dev = sigma / mean
        s = s + "{0:10}{1:^15.6E}{2:^15.6E}{3:^15.6E}{4:^15.6E}{5:^15.6E}{6:^15.6E}\n".format(iso, dat[iso][0], dat[iso][1], delta, mean, sigma, frac_dev)
    return s


#Open comparison file
cf = open(name + '_Compare.txt', 'w')

if type == "LWR":
    cf.write(header + "\n")

    p = subprocess.call("./LWR_PNL.py -p -c {0}".format(name), shell=True)

    pnl = 0.0
    with open(name + '_pnl.txt', 'r') as pnlf:
        pnl = [line.split()[0] for line in pnlf][0]
    cf.write("Calibrated to a Non-Leakage Probability of {0}\n\n".format(pnl))

    cf.write( SummaryTable("Input Isotopics:", -2) )
    cf.write("\n")
    cf.write( SummaryTable("Output Isotopics:", -1) )

    cf.write("\n\n\nHowever, when using a Non-Leakage Probability = 0.98\n")

    p = subprocess.call("./LWR_PNL.py -c {0}".format(name), shell=True)

    cf.write( SummaryTable("Input Isotopics:", -2) )
    cf.write("\n")
    cf.write( SummaryTable("Output Isotopics:", -1) )

    cf.write("\n\n\n\n")
elif type == "FR":
    cf.write(header + "\n")

    subprocess.call("./FR_PNL.py -c {0}".format(name), shell=True)

    pnl = 0.0
    with open(name + '_pnl.txt', 'r') as pnlf:
        pnl = [line.split()[0] for line in pnlf][0]
    cf.write("Calibrated to a Non-Leakage Probability of {0}\n".format(pnl))

    trucr = 0.0
    with open(name + '_trucr.txt', 'r') as trucrf:
        trucr = [line.split()[0] for line in trucrf][0]
    cf.write("Which gives a TRU_CR of {0}\n\n".format(trucr))

    cf.write( SummaryTable("Output Isotopics:", -1) )

    cf.write("\n\n\n\n")
else:
    print "Reactor type is wrong! Please pick either LWR or FR."

cf.close()
