import io
import os
import shutil
import logging

from bright.xsgen.plugins import Plugin
from bright.xsgen.testing import _run_tests
import bright.xsgen.envchar as envchar
from bright.xsgen.run.pbs import Pbs
from bright.xsgen.run.bash import Bash

from pyne.utils import failure, remove


class XDressPlugin(Plugin):
    """Core functionality for char."""

    # defaultrc = {"VERBOSE": False,
    #              "UI": False,
    #              "CWD": False,
    #              "CLEAN": False,
    #              "DEBUG": False,
    #              "TEST": False,
    #              "KILL_TRANSPORT": False,
    #              "PID": False,
    #              "FETCH_FILE": False,
    #              "MAKE_INPUT": False,
    #              "RUN_TRANSPORT": False,
    #              "RUN_BURNUP": False,
    #              "RUN_XS_GEN": False,
    #              "RUN_DELTAM": False,
    #              "RUN_ANALYSIS": False,
    # }

    def update_argparser(self, parser):
        # parser.add_argument("-v", "--verbose", action="store_true", dest="VERBOSE",
            # help="Gives extra info while running.")

        # parser.add_argument("--ui", action="store_true", dest="UI",
            # help="Launches the char ui.")

        # parser.add_argument("--cwd", action="store_true", dest="CWD",
        #     help="Run char in the current working directory.")

        # parser.add_argument("-c", "--clean", action="store_true", dest="CLEAN",
        #     help="Cleans the reactor directory of current files.")

        # parser.add_argument("--debug", action="store_true", dest="DEBUG",
        #     help="Turns on debug mode.")

        # parser.add_argument("-t", "--test", action="store_true", dest="TEST",
            # help="Tests an existing library for soundness.")

        # parser.add_argument("-k", "--kill", action="store_true", dest="KILL_TRANSPORT",
        #     help="Kills the current transport run. Sets -p.")

        # parser.add_argument("-p", "--pid", action="store_true", dest="PID",
        #     help="Finds the process identification number of a current transport run. Sets -d.")

        # parser.add_argument("-f", "--fetch", action="store_true", dest="FETCH_FILES",
        #     help="Fetches files from the remote server. Does not run transport, even if -r is set. Automatically sets -s.")

        # parser.add_argument("ABSPATH", help="A config file for xsgen.")

        # parser.add_argument("-i", "--input", action="store_true", dest="MAKE_INPUT",
            # help="Makes the transport calculation input deck.")

        # parser.add_argument("-r", "--run", action="store_true", dest="RUN_TRANSPORT",
            # help="Run the transport calculation.")

        # parser.add_argument("-b", "--burnup",  action="store_true", dest="RUN_BURNUP",
        #     help="Run the burnup calculation.")

        # parser.add_argument("-x", "--xs",  action="store_true", dest="RUN_XS_GEN", 
        #     help="Run the cross-section generation calculation.")

        # parser.add_argument("-m", "--delta-mass",  action="store_true", dest="RUN_DELTAM", 
        #     help="Run the initial nuclide sensitivity calculation.")

        # parser.add_argument("-a", "--analyze", action="store_true", dest="RUN_ANALYSIS", 
        #     help="Run analysis on database.")


    def execute(self, rc):

        # run_switch = {'': Bash,
                      # 'BASH': Bash,
                      # 'bash': Bash,
                      # 'PBS': Pbs,
                      # 'pbs': Pbs,
                      # 'Torque': Pbs,
                      # 'torque': Pbs,
                      # }

        # from n_code_serpent import NCodeSerpent
        # from mock import Mock
        # NCodeSerpent = Mock()

        # n_code_switch = {'': NCodeSerpent,
        #                  'sss': NCodeSerpent,
        #                  'Serpent': NCodeSerpent,
        #                  'serpent': NCodeSerpent,
        #                  }

        # if rc.UI:
        #     # Test to see if ui library is installed
        #     try:
        #         from bright.xsgen.ui import app
        #     except ImportError:
        #         print(failure("Please install the Enthought Tool Suite (ETS) for CHAR UI."))
        #         raise SystemExit

        #     # Open UI
        #     application = app.Application()
        #     #application.rx_h5_path = "/home/scopatz/MultiGroupPaper/DataXS/lwr/lwr.h5"
        #     application.configure_traits()

        #     # Clean-up UI
        #     if application.rx_h5 is not None:
        #         application.rx_h5.close()

        #     raise SystemExit

        # # Load the CHAR definition file into the RC
        # execfile(rc.ABSPATH, {}, rc._dict)

        # # update defchar adding more useful values.
        # rc = envchar.update_env(rc)

        # if rc.TEST:
            # _run_tests(rc.ABSPATH)
            # raise SystemExit

        # if rc.KILL_TRANSPORT:
        #     rc.PID = True                      #Ensures that the PID is found in order that it mak be killed.

        # if rc.PID:
        #     rc.RUN_XS_GEN = False
        #     rc.RUN_DELTAM = False
        #     rc.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.

        # if rc.FETCH_FILES:
        #     rc.RUN_XS_GEN = False
        #     rc.RUN_DELTAM = False
        #     rc.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.
        #     rc.LOCAL = False                   #Ensures that ssh package is loaded.

        # if not rc.CWD:
        #     if rc.CLEAN:
        #         remove(rc.reactor)
        #     if rc.reactor not in os.listdir('.'):
        #         os.mkdir(rc.reactor)

        #     os.chdir(rc.reactor)
        #     shutil.copyfile(rc.ABSPATH, 'defchar.py')

        # # Start up logger
        # logger = logging.getLogger('char')
        # hdlr = logging.FileHandler('char.log')
        # formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        # hdlr.setFormatter(formatter)
        # logger.addHandler(hdlr)
        # logger.setLevel(logging.INFO)
        # rc.logger = logger

        # # Set the neutronics code
        # n_coder = n_code_switch[rc.transporter]
        # n_code = n_coder(rc)
        # rc.n_code = n_code

        # # Get the run controller
        # runner = run_switch[rc.scheduler]
        # runchar = runner(n_code, rc)

        # if (rc.MAKE_INPUT) and (not rc.FETCH_FILES) and (not rc.PID):
            # runchar.init_h5()

        # # Check a bunch of run conditions
        # if rc.RUN_TRANSPORT:
        #     # Run Transport code
        #     runchar.make_run_script()

        #     if rc.LOCAL:
        #         runchar.run_locally()
        #     else:
        #         runchar.run_remotely()

        # elif rc.RUN_ANALYSIS:
        #     n_code.analyze_deltam()

        # elif rc.RUN_BURNUP or rc.RUN_XS_GEN or rc.RUN_DELTAM:
        #     # Make tranumatrion libraries by executing the as a separate step from 
        #     # the cross-section generation
        #     if rc.RUN_BURNUP:
        #         runchar.burnup(idx)

        #     # Make Cross-sections as a separate step from the burnup calculation
        #     if rc.RUN_XS_GEN:
        #         runchar.xs_gen(idx, nucs)

        #     # Run initial nuclide sensitivity calculation
        #     if rc.RUN_DELTAM:
        #         n_code.run_deltam_pert(idx, ihm_isos, sidx)

        # elif rc.FETCH_FILES:
        #     # Fetches files from remote server
        #     runchar.fetch()
        # elif rc.PID:
        #     # Finds and prints the PID of CHAR
        #     runchar.pid()
        # elif rc.KILL_TRANSPORT:
        #     # Finds and kills CHAR
        #     runchar.kill()

        # Clean up
        # if not rc.CWD:
            # os.chdir('..')
