from plugins import Plugin
import os

# from utils import RunControl, NotSpecified, DEFAULT_RC_FILE, DEFAULT_PLUGINS

class XDressPlugin(Plugin):
    """Core functionality for char."""

    def update_argparser(self, parser):
        parser.add_option("--cwd", action="store_true", dest="CWD", default=False, 
            help="Run char in the current working directory.")

        parser.add_option("-c", "--clean", action="store_true", dest="CLEAN", 
            help="Cleans the reactor direactory of current files.")

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
                os.remove(rc.reactor)

            if rc.reactor not in os.listdir('.'):
                os.mkdir(rc.reactor)

            os.chdir(rc.reactor)
            shutil.copyfile(rc.absolute_path, 'defchar.py')
