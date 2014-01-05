#! /usr/bin/env python2.7

############################
#### Standard Libraries ####
############################
from __future__ import print_function
import os
import shutil
import logging

from optparse import OptionParser

###########################
### Extension Libraries ###
###########################

##########################
#### Custom Libraries ####
##########################

from pyne import nucname
name_zz = nucname.name_zz
zz_name = nucname.zz_name

#################
### CHAR Libs ###
#################
from bright.xsgen import envchar

from bright.xsgen.run.pbs import Pbs
from bright.xsgen.run.bash import Bash

from pyne.utils import failure, remove

from bright.xsgen.testing import _run_tests

run_switch = {'': Bash, 
              'BASH': Bash,
              'bash': Bash,
              'PBS': Pbs,
              'pbs': Pbs,
              'Torque': Pbs,
              'torque': Pbs,
              }

# from n_code_serpent import NCodeSerpent
from mock import Moc
NCodeSerpent = Mock()

n_code_switch = {'': NCodeSerpent, 
                 'sss': NCodeSerpent, 
                 'Serpent': NCodeSerpent, 
                 'serpent': NCodeSerpent, 
                 }


def parse_slice(s, size):
    """Parses a string into a list of indices."""
    l = [int(i) for i in s.split(':') if 0 < len(i)]

    # Handle negative indices
    for n in range(len(l)):
        if l[n] < 0:
            l[n] = size + l[n]

    if len(l) == 0:
        l = [0, size]

    if len(i) == 1:
        l.append(l[0] + 1)

    return l


def parse_nucs(s):
    """Parses a string into a set of nuclides."""
    nset = set()
    nucs = s.split(',')

    for nuc in nucs:
        if len(nuc) == 0:
            continue
        elif '-' in nuc:
            nsplit = nuc.split()
            nlower = nucname.zzaaam(nsplit[0])
            nupper = nucname.zzaaam(nsplit[1])
            if 0 == nupper%10000:
                nupper += 10000
            else:
                nupper += 1
            tmpset = set(range(nlower, nupper))
        else:
            n = nucname.zzaaam(nuc)
            if 0 == n%10000:
                nrange = range(n, n + 10000)
            else:
                nrange = [n]
            tmpset = set(nrange)

        # Add the union 
        nset = (nset | tmpset)

    return nset


def main():
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

    parser.add_option("-a", "--analyze", action="store_true", dest="RUN_ANALYSIS", 
        default=False, help="Run analysis on database.")

    parser.add_option("-b", "--burnup",  action="store_true", dest="RUN_BURNUP", 
        default=False, help="Run the burnup calculation.")

    parser.add_option("-x", "--xs",  action="store_true", dest="RUN_XS_GEN", 
        default=False, help="Run the cross-section generation calculation.")

    parser.add_option("-m", "--delta-mass",  action="store_true", dest="RUN_DELTAM", 
        default=False, help="Run the initial nuclide sensitivity calculation.")

    parser.add_option("-N", action="store", dest="NPERT", 
        default='', help="Pertubation indices to calculate, in Python slice syntax.")

    parser.add_option("-I",  action="store", dest="ISOS", 
        default='', help="Nuclides to calculate.")

    parser.add_option("-S",  action="store", dest="NSENS", 
        default='', help="Sensitivity indices to calculate, in Python slice syntax.")

    parser.add_option("-c", "--clean", action="store_true", dest="CLEAN", 
        help="Cleans the reactor direactory of current files.")

    parser.add_option("-C", "--cache", action="store_true", dest="CACHE", 
        default=False, help="Uses the current files in the reactor direactory.")

    parser.add_option("-l", "--local", action="store_true", dest="LOCAL", 
        default=True, help="Run or Fetch files locally.")

    parser.add_option("-s", "--server", action="store_false", dest="LOCAL", 
        help="Run or Fetch files from a remote server.")

    parser.add_option("-f", "--fetch", action="store_true", dest="FETCH_FILES", default=False, 
        help="Fetches files from the remote server. Does not run transport, even if -r is set. Automatically sets -s.")

    parser.add_option("-p", "--pid", action="store_true", dest="PID", default=False, 
        help="Finds the process identification number of a current transport run. Sets -d.")

    parser.add_option("-k", "--kill", action="store_true", dest="KILL_TRANSPORT", 
        default=False, help="Kills the current transport run. Sets -p.")

    parser.add_option("--cwd", action="store_true", dest="CWD", default=False, 
        help="Run char in the current working directory.")

    parser.add_option("--ui", action="store_true", dest="UI", default=False, 
        help="Launches the char ui.")

    parser.add_option("-t", "--test", action="store_true", dest="TEST", 
        default=False, help="Tests an existing library for soundness.")

    (options, args) = parser.parse_args()

    # Try launching ui before anything else
    if options.UI:
        # Test to see if ui library is installed
        try:
            from .ui import app
        except ImportError:
            print(failure("Please install the Enthought Tool Suite (ETS) for CHAR UI."))
            raise SystemExit

        # Open UI
        application = app.Application()
        #application.rx_h5_path = "/home/scopatz/MultiGroupPaper/DataXS/lwr/lwr.h5"
        application.configure_traits()

        # Clean-up UI
        if application.rx_h5 is not None:
            application.rx_h5.close()

        raise SystemExit

    # Make sure we have a configureation file before proceeding
    if len(args) == 0:
        print(failure("Please specify a file for char."))
        raise SystemExit

    absolute_path = os.path.abspath(args[0])

    # Run tests on a db file
    if options.TEST:
        _run_tests(absolute_path)
        raise SystemExit

    # Load the CHAR definition file into its own env namespace
    env = {}
    execfile(absolute_path, {}, env)


    # Add command line arguments to env
    env['options'] = options
    env['args'] = args

    # Update defchar adding more useful values.
    env = envchar.update_env(env)

    #intial command-line options protocol.
    if options.KILL_TRANSPORT:
        options.PID = True                      #Ensures that the PID is found in order that it mak be killed.

    if options.PID:
        options.RUN_XS_GEN = False
        options.RUN_DELTAM = False
        options.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.

    if options.FETCH_FILES:
        options.RUN_XS_GEN = False
        options.RUN_DELTAM = False
        options.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.
        options.LOCAL = False                   #Ensures that ssh package is loaded.

    ################
    #### Script ####
    ################

    # Prep work
    if not options.CWD:
        if env['options'].CLEAN:
            remove(env['reactor'])

        if env['reactor'] not in os.listdir('.'):
            os.mkdir(env['reactor'])

        os.chdir(env['reactor'])
        shutil.copyfile(absolute_path, 'defchar.py')

    # Start up logger
    logger = logging.getLogger('char')
    hdlr = logging.FileHandler('char.log')
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)
    env['logger'] = logger

    # Set the neutronics code
    n_coder = n_code_switch[env['transporter']]
    n_code = n_coder(env)
    env['n_code'] = n_code

    # Get the run controller
    runner = run_switch[env['scheduler']]
    runchar = runner(n_code, env)

    # Get the inputs indeces
    idx = parse_slice(options.NPERT, n_code.nperturbations)

    nucs = parse_nucs(options.ISOS)
    if len(nucs) == 0:
        nucs = set(env['core_transmute'])
    else:
        nucs = (nucs & set(env['core_transmute']))
    ihm_nucs = (nucs & set(env['ihm_mat'].comp.keys()))

    if 'deltam' in env:
        sidx = parse_slice(options.NSENS, len(env['deltam']))

    # Make the input file unless otherwise specified.
    if (options.MAKE_INPUT) and (not options.FETCH_FILES) and (not options.PID):
        runchar.init_h5()

    # Check a bunch of run conditions
    if options.RUN_TRANSPORT:
        # Run Transport code
        runchar.make_run_script()

        if options.LOCAL:
            runchar.run_locally()
        else:
            runchar.run_remotely()

    elif options.RUN_ANALYSIS:
        n_code.analyze_deltam()

    elif options.RUN_BURNUP or options.RUN_XS_GEN or options.RUN_DELTAM:
        # Make tranumatrion libraries by executing the as a separate step from 
        # the cross-section generation
        if options.RUN_BURNUP:
            runchar.burnup(idx)

        # Make Cross-sections as a separate step from the burnup calculation
        if options.RUN_XS_GEN:
            runchar.xs_gen(idx, nucs)

        # Run initial nuclide sensitivity calculation
        if options.RUN_DELTAM:
            n_code.run_deltam_pert(idx, ihm_isos, sidx)

    elif options.FETCH_FILES:
        # Fetches files from remote server
        runchar.fetch()
    elif options.PID:
        # Finds and prints the PID of CHAR
        runchar.pid()
    elif options.KILL_TRANSPORT:
        # Finds and kills CHAR
        runchar.kill()

    # Clean up
    if not options.CWD:
        os.chdir('..')


# Run CHAR
if __name__ == '__main__':
    main()
