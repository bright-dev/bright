import argparse
from bright.xsgen.plugins import Plugins
from bright.xsgen.utils import NotSpecified, RunControl, exec_file
import os

try:
    import argcomplete
except ImportError:
    argcomplete = None

DEFAULT_RC_FILE="rc"
DEFAULT_PLUGINS=["base"]

def main():
    preparser = argparse.ArgumentParser()
    preparser.add_argument("infile")
    preparser.add_argument("--plugins", default=NotSpecified, nargs="+")
    preparser.add_argument("--rc", default=NotSpecified, nargs="+")
    prens = preparser.parse_known_args()[0]
    prens.abspath = os.path.abspath(prens.infile)
    predefaultrc = RunControl(rc=DEFAULT_RC_FILE, plugins=DEFAULT_PLUGINS)
    prerc = RunControl()
    prerc._update(predefaultrc)
    prerc.rc = prens.rc
    rcdict = {}
    if os.path.isfile(prerc.rc):
        exec_file(prerc.rc, rcdict, rcdict)
        prerc.rc = rcdict['rc'] if 'rc' in rcdict else NotSpecified
        prerc.plugins = rcdict['plugins'] if 'plugins' in rcdict else NotSpecified
    prerc._update([(k, v) for k, v in prens.__dict__.items()])

    plugins = Plugins(prerc.plugins)
    parser = plugins.build_cli()
    if argcomplete is not None and prerc.bash_completion:
        argcomplete.autocomplete(parser)
    ns = parser.parse_args()
    rc = plugins.merge_rcs()
    rc._update(rcdict)
    rc._update([(k, v) for k, v in ns.__dict__.items()])
    plugins.setup()
    plugins.execute()
    plugins.teardown()


if __name__ == "__main__":
    main()
