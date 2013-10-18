from plugins import Plugin
# from utils import RunControl, NotSpecified, DEFAULT_RC_FILE, DEFAULT_PLUGINS

class XDressPlugin(Plugin):
    """Core functionality for char."""
    # defaultrc = RunControl(
    #     rc=DEFAULT_RC_FILE,
    #     plugins=DEFAULT_PLUGINS,
    #     debug=False,
    #     ts=TypeSystem(),
    #     verbose=False,
    #     dumpdesc=False,
    #     package=NotSpecified,
    #     packagedir=NotSpecified,
    #     sourcedir='src',
    #     builddir='build',
    #     bash_completion=True,
    #     )

    # Sweet hack because ts.update() returns None
    rcupdaters = {'ts': (lambda old, new: old.update(new) or old)}

    rcdocs = {
        'rc': "Path to run control file",
        'plugins': "Plugins to include",
        'debug': 'Build in debugging mode',
        'verbose': "Print more output.",
        'ui': "Use the GUI."
        }

    def update_argparser(self, parser):
        parser.add_argument('--rc', help=self.rcdocs['rc'])
        parser.add_argument('--plugins', nargs="+", help=self.rcdocs["plugins"])
        parser.add_argument('--debug', action='store_true',
                            help=self.rcdocs["debug"])
        parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
                            help=self.rcdocs["verbose"])


    def report_debug(self, rc):
        msg = 'Current descripton cache contents:\n\n{0}\n\n'
        msg = msg.format(str(rc._cache))
        msg += nyansep + "\n\n"
        msg += "Current type system contents:\n\n" + str(rc.ts) + "\n\n"
        return msg

