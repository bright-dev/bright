"""The base plugin for xsgen.

Base Plugin API
===============
"""
from __future__ import print_function
import os
import sys
from warnings import warn

from bright.xsgen.utils import RunControl, NotSpecified, writenewonly, \
    DEFAULT_RC_FILE, DEFAULT_PLUGINS, nyansep, indent
from bright.xsgen.plugins import Plugin
from bright.xsgen.version import report_versions

if sys.version_info[0] >= 3:
    basestring = str

class XSGenPlugin(Plugin):
    """This class provides base functionality for xdress itself."""

    defaultrc = RunControl(
        rc=DEFAULT_RC_FILE,
        plugins=DEFAULT_PLUGINS,
        debug=False,
        verbose=False,
        version=False,
        bash_completion=True,
        )

    rcdocs = {
        'rc': "Path to run control file",
        'plugins': "Plugins to include",
        'debug': 'Build in debugging mode', 
        'verbose': "Print more output.",
        'version': "Print version information.",
        'bash_completion': ("Flag for enabling / disabling BASH completion. "
                            "This is only relevant when using argcomplete."),
        }

    def update_argparser(self, parser):
        parser.add_argument('--rc', help=self.rcdocs['rc'])
        parser.add_argument('--plugins', nargs="+", help=self.rcdocs["plugins"])
        parser.add_argument('--debug', action='store_true', 
                            help=self.rcdocs["debug"])
        parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
                            help=self.rcdocs["verbose"])
        parser.add_argument('--version', action='store_true', dest='version',
                            help=self.rcdocs["version"])
        parser.add_argument('--bash-completion', action='store_true',
                            help="enable bash completion", dest="bash_completion")

    def setup(self, rc):
        if rc.version:
            print(report_versions())
            sys.exit()

    def report_debug(self, rc):
        msg = 'Version Information:\n\n{0}\n\n'
        msg = msg.format(indent(report_versions()))
        msg += nyansep + "\n\n"
        return msg

