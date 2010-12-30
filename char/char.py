#!/usr/bin/env python

############################
#### Standard Libraries ####
############################
from __future__ import print_function
import os
import sys
import time
import shutil
import logging
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

from metasci.colortext import failure

#################
### CHAR Libs ###
#################
defchar = None
import glbchar

import graphchar
import runchar
from n_code_mcnp    import NCodeMCNP
from n_code_origen  import NCodeORIGEN
from n_code_serpent import NCodeSerpent

def main():
    global defchar

    ###########################
    ### Command Line Parser ###
    ###########################
    usage = "usage: %prog [options] confchar"
    parser = OptionParser(usage)

    parser.add_option("-v", "--verbose", action="store_true", dest="VERBOSE", 
        default=False, help="Gives extra info while running.")

    parser.add_option("-i", "--input", action="store_true", dest="MAKE_INPUT", 
        help="Makes the transport calculation input deck.")

    parser.add_option("-r", "--run", action="store_true", dest="RUN_TRANSPORT", 
        default=False, help="Run the transport calculation.")

    parser.add_option("-d", "--dry-run", action="store_false", dest="RUN_TRANSPORT", 
        help="Dry Run. Do NOT run the transport calculation.")

    parser.add_option("-b", "--burnup",  action="store_true", dest="RUN_BURNUP", 
        default=False, help="Run the burnup calculation.")

    parser.add_option("-x", "--xs",  action="store_true", dest="RUN_XS_GEN", 
        default=False, help="Run the cross-section generation calculation.")

    parser.add_option("-c", "--clean", action="store_true", dest="CLEAN", 
        help="Cleans the reactor direactory of current files.")

    parser.add_option("-l", "--local", action="store_true", dest="LOCAL", 
        default=True, help="Run or Fetch files locally.")

    parser.add_option("-s", "--server", action="store_false", dest="LOCAL", 
        help="Run or Fetch files from a remote server.")

    parser.add_option("-f", "--fetch", action="store_true", dest="FETCH_FILES", default=False, 
        help="Fetches files from the remote server. Does not run transport, even if -r is set. Automatically sets -s.")

    parser.add_option("-p", "--pid", action="store_true", dest="PID", default=False, 
        help="Finds the process identification number of a current transport run. Sets -d.")

    parser.add_option("-k", "--kill", action="store_true", dest="KILL_TRANSPORT", 
        default=False, help="Kills the current transport run. Sets -P.")

    parser.add_option("--cwd", action="store_true", dest="CWD", default=False, 
        help="Run char in the current working directory.")

    (options, args) = parser.parse_args()

    # Make sure we have a configureation file before proceeding
    if len(args) == 0:
        print(failure("Please specify a configuration file for CHAR."))
        raise SystemExit

    # Load the CHAR definition file early, into its own namespace
    absolute_path = os.path.abspath(args[0])
    dir, file = os.path.split(absolute_path)
    sys.path.append(dir)
    mod_name = file.rpartition('.')[0]
    defchar = __import__(mod_name)

    # Add command line arguments to defchar
    defchar.options = options
    defchar.args = args

    # Update defchar adding more useful values.
    defchar = glbchar.defchar_update(defchar)

    #intial command-line options protocol.
    if options.KILL_TRANSPORT:
        options.PID = True                      #Ensures that the PID is found in order that it mak be killed.

    if options.PID:
        options.RUN_XS_GEN = False
        options.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.

    if options.FETCH_FILES:
        options.RUN_XS_GEN = False
        options.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.
        options.LOCAL = False                   #Ensures that ssh package is loaded.

    ################
    #### Script ####
    ################

    # Prep work
    if not options.CWD:
        if defchar.options.CLEAN:
            metasci.safe_remove(defchar.reactor, True)

        if defchar.reactor not in os.listdir('.'):
            os.mkdir(defchar.reactor)

        os.chdir(defchar.reactor)
        shutil.copyfile(absolute_path, 'defchar.py')

    # Start up logger
    logger = logging.getLogger('char')
    hdlr = logging.FileHandler('char.log')
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)
    defchar.logger = logger

    # Set the transport code type
    try:
        transport_code = defchar.transport_code.lower()
        defchar.transport_code = transport_code
    except:
        print(failure("Transport code type not set!\n"
                      "  Use either the '-w' command line option, or\n"
                      "  use the 'TransportCode' flag in defchar.py.\n"
                      "Currently, 'MCNP' and 'Serpent' are accpeted values.\n"
                      "\n"
                      "Note: you must have the appropriate neutronics code\n"
                      "installed on your machine for this to work."))
        raise SystemExit

    if ('serpent' in transport_code) and ('mcnp' not in transport_code):
        n_transporter = NCodeSerpent()
    elif ('mcnp' in transport_code) and ('serpent' not in transport_code):
        n_transporter = NCodeMCNP()
    else:
        print(failure("The transport code given is not valid: {0:yellow}", transport_code))
        print(failure("Currently, 'MCNP' and 'Serpent' are accpeted values."))
        raise SystemExit

    defchar.n_transporter = n_transporter

    # Make the input file unless otherwise specified.
    if (options.MAKE_INPUT) and (not options.FETCH_FILES) and (not options.PID):
        n_transporter.make_input()

    # Check a bunch of run conditions
    if options.RUN_TRANSPORT:
        # Run Transport code
        runchar.make_run_script(n_transporter)

        if options.LOCAL:
            runchar.run_transport_local()
        else:
            runchar.run_transport_remote()

    elif options.RUN_BURNUP or options.RUN_XS_GEN:
        # Make tranumatrion libraries by executing the as a separate step from 
        # the cross-section generation
        if options.RUN_BURNUP:
            n_transporter.run_burnup()

        # Make Cross-sections as a separate step from the burnup calculation
        if options.RUN_XS_GEN:
            n_transporter.run_xs_gen()

    elif options.FETCH_FILES:
        #Fetches files from remote server
        runchar.fetch_remote_files()
    elif options.PID:
        #Finds (and kills?) the Transport Run Process
        if options.LOCAL:
            runchar.find_pid_local()
        else:
            runchar.find_pid_remote()

    #Clean up
    if not options.CWD:
        os.chdir('..')


#Run CHAR
if __name__ == '__main__':
    main()
