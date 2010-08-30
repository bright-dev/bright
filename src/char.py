#!/usr/bin/env python

############################
#### Standard Libraries ####
############################
from __future__ import print_function
import os
import time
import shutil
import subprocess

from optparse import OptionParser

###########################
### Extension Libraries ###
###########################
import tables as tb

##########################
#### Custom Libraries ####
##########################
import isoname
import metasci
import metasci.nuke as msn
import metasci.graph as msg

from metasci.colortext import *

#################
### CHAR Libs ###
#################
from glbchar import *

import parsechar
import graphchar
import runchar

def main():
    ###########################
    ### Command Line Parser ###
    ###########################
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option("-q", "--quiet", action="store_true", dest="Quiet", 
        default=False, help="Supresses stdout.")

    parser.add_option("-v", "--verbose", action="store_true", dest="Verbose", 
        default=False, help="Gives extra info while running.")

    parser.add_option("-i", "--input", action="store_true", dest="MakeInput", 
        help="Makes the MCNP input deck.")

    parser.add_option("-r", "--run", action="store_true", dest="RunTransport", 
        default=False, help="Run the transport calculation.")

    parser.add_option("-d", "--dry-run", action="store_false", dest="RunTransport", 
        help="Dry Run. Do NOT run the transport calculation.")

    parser.add_option("-w", "--with", dest="RunWith", default="NONE", metavar="PROG", 
        help="Dictates what PROG to run transport calculation with. PROG = [MCNP | MPI | PBS].")

    parser.add_option("-O", "--ORIGEN", action="store_true", dest="RunORIGEN", 
        default=False, help="Run ORIGEN Burnup calculations.")

    parser.add_option("-l", "--local", action="store_true", dest="Local", 
        default=True, help="Run or Fetch files locally.")

    parser.add_option("-s", "--server", action="store_false", dest="Local", 
        help="Run or Fetch files from a remote server.")

    parser.add_option("-f", "--fetch", action="store_true", dest="FetchFiles", default=False, 
        help="Fetches files from the remote server. Does not run transport, even if -r is set. Automatically sets -s.")

    parser.add_option("-p", "--parse", action="store_true", dest="ParseData", 
        default=False, help="Parses output and stores it an HDF5 library.")

    parser.add_option("-t", "--text", action="store_true", dest="MakeText", 
        default=False, help="Parses output and makes text-based libraries.")

    parser.add_option("-9", "--tape9", action="store_true", dest="MakeTape9", 
        default=False, help="Parses output and makes a ORIGEN tape9 libraries for each burn step.")

    parser.add_option("-g", "--graph", action="store_true", dest="MakeGraphs", 
        default=False, help="Graphs data from libraries. Sets -p.")

    parser.add_option("-G", "--graph-ORIGEN", action="store_true", dest="MakeOrigenGraphs", 
        default=False, help="Graphs ORIGEN data only.")

    parser.add_option("-P", "--pid", action="store_true", dest="PID", default=False, 
        help="Finds the process identification number of a current transport run. Sets -d.")

    parser.add_option("-k", "--kill", action="store_true", dest="KillTransport", 
        default=False, help="Kills the current transport run. Sets -P.")

    parser.add_option("--no-burn", action="store_true", dest="NoBurnBool", default=False, 
        help="Removes the burn card from the MCNP deck.")

    parser.add_option("--no-pert", action="store_true", dest="NoPertBool", default=False, 
        help="Removes the perturbation cards from the MCNP deck.")

    (options, args) = parser.parse_args()

    #intial command-line options protocol.
    Quiet = options.Quiet
    options.RunWith = options.RunWith.upper()   #Ensures that RunWith is uppercase (nice for flags).
    if options.KillTransport:
        options.PID = True                      #Ensures that the PID is found in order that it mak be killed.
    if options.PID:
        options.RunTransport = False            #Ensures that transport calculation is not initiated while fetching files.
    if options.FetchFiles:
        options.RunTransport = False            #Ensures that transport calculation is not initiated while fetching files.
        options.Local = False                   #Ensures that ssh package is loaded.
    if options.RunORIGEN:
        options.MakeTape9 = True
    if options.RunWith == "NONE":               #Sets defualt run behaviour
        if options.Local:
            options.RunWith = "MCNP"
        else:
            options.RunWith = "PBS"

    ################
    #### Script ####
    ################

    #Prep work
    os.chdir(reactor)

    #Make the input file unless otherwise specified.
    if (options.MakeInput) and (not options.FetchFiles) and (not options.PID):
        runchar.Make_MCNP_Input(options.NoBurnBool, options.NoPertBool)

    #Run MCNP
    if options.RunTransport:
        metasci.SafeRemove(reactor + ".o")
        metasci.SafeRemove(reactor + ".m")
        metasci.SafeRemove(reactor + ".r")
        metasci.SafeRemove(reactor + ".s")

        runchar.Make_Run_Script(options.RunWith, options.Local)
        if options.Local:
            runchar.Run_Transport_Local(options.RunWith)
        else:
            runchar.Run_Transport_Remote(options.RunWith)

    #Fetches files from remote server
    if options.FetchFiles:
        runchar.Fetch_Remote_Files()

    #Finds (and kills?) the Transport Run Process
    if options.PID:
        if options.Local:
            runchar.Find_PID_Local(options.KillTransport)
        else:
            runchar.Find_PID_Remote(options.RunWith, options.KillTransport)

    #Parse MCNPX Output & Make HDF5 Data Library
    if options.ParseData:
        parsechar.Write_HDF5_Lib_MCNP()

        #Parse MCNPX Output & Make Text-based Data Libraries
        if options.MakeText:
            parsechar.Write_TXT_Lib_MCNP()

    #Parse MCNPX Output & Make ORIGEN TAPE9 Libraries
    if options.MakeTape9:
        parsechar.Write_ORIGEN_Libs()

    #Run ORIGEN IRF Burnup Calculation
    if options.RunORIGEN:
        BU, k, Pro, Des, Tij = runchar.Run_ORIGEN()
        parsechar.Write_HDF5_Lib_ORIGEN( BU, k, Pro, Des, Tij )

        if options.MakeText:
            parsechar.Write_TXT_Lib_ORIGEN( BU, k, Pro, Des, Tij )

    #Make Graphs, for fun and learning
    if options.MakeGraphs:
        graphchar.Make_Figs_MCNP(options.Quiet)

        if options.RunORIGEN:
            graphchar.Make_Figs_ORIGEN(options.Quiet)

    #Make only the ORIGEN figures
    if options.MakeOrigenGraphs and not (options.MakeGraphs and options.RunORIGEN):
        graphchar.Make_Figs_ORIGEN(options.Quiet)


    #Clean up
    metasci.SafeRemove(reactor + ".m")
    metasci.SafeRemove(reactor + ".r")
    metasci.SafeRemove(reactor + ".s")
    os.chdir('..')


#Run CHAR
if __name__ == '__main__':
    main()
