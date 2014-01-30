import logging
import os
import shutil

from bright.xsgen.plugins import Plugin
import bright.xsgen.envchar as envchar
from bright.xsgen.run.pbs import Pbs
from bright.xsgen.run.bash import Bash

class XDressPlugin(Plugin):

    defaultrc = {"VERBOSE": False,
                 "UI": False,
                 "CWD": False,
                 "CLEAN": False,
                 "DEBUG": False,
                 "TEST": False,
                 "KILL_TRANSPORT": False,
                 "PID": False,
                 "FETCH_FILE": False,
                 "MAKE_INPUT": False,
                 "RUN_TRANSPORT": False,
                 "RUN_BURNUP": False,
                 "RUN_XS_GEN": False,
                 "RUN_DELTAM": False,
                 "RUN_ANALYSIS": False,
    }

    def update_argparser(self, parser):
        parser.add_argument("--ui", action="store_true", dest="UI",
            help="Launches the char ui.")

        parser.add_argument("--cwd", action="store_true", dest="CWD",
            help="Run char in the current working directory.")

        parser.add_argument("-c", "--clean", action="store_true", dest="CLEAN",
            help="Cleans the reactor directory of current files.")

        parser.add_argument("-v", "--verbose", action="store_true", dest="VERBOSE",
            help="Gives extra info while running.")

        parser.add_argument("--debug", action="store_true", dest="DEBUG",
            help="Turns on debug mode.")

        parser.add_argument("ABSPATH", help="A config file for xsgen.")

        parser.add_argument("--rc", help="Alternative rc file")

    def setup(self, rc):
        run_switch = {'': Bash,
                      'BASH': Bash,
                      'bash': Bash,
                      'PBS': Pbs,
                      'pbs': Pbs,
                      'Torque': Pbs,
                      'torque': Pbs,
                      }

        # Get the run controller
        runner = run_switch[rc.scheduler]
        rc.runchar = runner(n_code, rc)

        from mock import Mock
        NCodeSerpent = Mock()

        # Set the neutronics code
        n_code_switch = {'': NCodeSerpent,
                         'sss': NCodeSerpent,
                         'Serpent': NCodeSerpent,
                         'serpent': NCodeSerpent,
                         }
        n_coder = n_code_switch[rc.transporter]
        n_code = n_coder(rc)
        rc.n_code = n_code

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

        # Load the CHAR definition file into the RC
        execfile(rc.ABSPATH, {}, rc._dict)

        # update defchar adding more useful values.
        rc = envchar.update_env(rc)

        if not rc.CWD:
            if rc.CLEAN:
                remove(rc.reactor)
            if rc.reactor not in os.listdir('.'):
                os.mkdir(rc.reactor)
            os.chdir(rc.reactor)
            shutil.copyfile(rc.ABSPATH, 'defchar.py')

        # Start up logger
        logger = logging.getLogger('char')
        hdlr = logging.FileHandler('char.log')
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr)
        logger.setLevel(logging.INFO)
        rc.logger = logger
