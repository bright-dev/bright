from bright.xsgen.plugins import Plugin
from bright.xsgen.testing import _run_tests

class XDressPlugin(Plugin):
    def update_argparser(self, parser):
        parser.add_argument("-t", "--test", action="store_true", dest="TEST",
            help="Tests an existing library for soundness.")

    def execute(self, rc):
        if rc.TEST:
            _run_tests(rc.ABSPATH)
            raise SystemExit
