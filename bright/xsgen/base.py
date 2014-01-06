import io
import os
import shutil

from bright.xsgen.plugins import Plugin
from bright.xsgen.testing import _run_tests
import bright.xsgen.envchar as envchar

from pyne.utils import failure, remove


class XDressPlugin(Plugin):
    """Core functionality for char."""

    defaultrc = {"UI": False,
                 "CWD": False,
                 "CLEAN": False,
                 "DEBUG": False,
                 "TEST": False,
                 "KILL_TRANSPORT": False,
                 "PID": False,
                 "FETCH_FILE": False,
    }
    def update_argparser(self, parser):
        parser.add_argument("--ui", action="store_true", dest="UI",
            help="Launches the char ui.")

        parser.add_argument("--cwd", action="store_true", dest="CWD",
            help="Run char in the current working directory.")

        parser.add_argument("-c", "--clean", action="store_true", dest="CLEAN",
            help="Cleans the reactor directory of current files.")

        parser.add_argument("--debug", action="store_true", dest="DEBUG",
            help="Turns on debug mode.")

        parser.add_argument("-t", "--test", action="store_true", dest="TEST",
            help="Tests an existing library for soundness.")

        parser.add_argument("-k", "--kill", action="store_true", dest="KILL_TRANSPORT",
            help="Kills the current transport run. Sets -p.")

        parser.add_argument("-p", "--pid", action="store_true", dest="PID",
            help="Finds the process identification number of a current transport run. Sets -d.")

        parser.add_argument("-f", "--fetch", action="store_true", dest="FETCH_FILES",
            help="Fetches files from the remote server. Does not run transport, even if -r is set. Automatically sets -s.")

        parser.add_argument("ABSPATH", help="A config file for xsgen.")

    def execute(self, rc):
        if rc.UI:
            # Test to see if ui library is installed
            try:
                from bright.xsgen.ui import app
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

        if not rc.CWD:
            if rc.CLEAN:
                remove(rc.reactor)

            if rc.reactor not in os.listdir('.'):
                os.mkdir(rc.reactor)

            os.chdir(rc.reactor)
            shutil.copyfile(rc.absolute_path, 'defchar.py')

        if rc.TEST:
            _run_tests(rc.ABSPATH)
            raise SystemExit

        # Load the CHAR definition file into the RC
        execfile(rc.ABSPATH, {}, rc)

        # update defchar adding more useful values.
        rc = envchar.update_env(rc)

        if rc.KILL_TRANSPORT:
            rc.PID = True                      #Ensures that the PID is found in order that it mak be killed.

        if rc.PID:
            rc.RUN_XS_GEN = False
            rc.RUN_DELTAM = False
            rc.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.

        if rc.FETCH_FILES:
            rc.RUN_XS_GEN = False
            rc.RUN_DELTAM = False
            rc.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.
            rc.LOCAL = False                   #Ensures that ssh package is loaded.
