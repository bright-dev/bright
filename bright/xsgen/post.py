from bright.xsgen.plugins import Plugin

class XDressPlugin(Plugin):
    def execute(self, rc):
        # Clean up
        if not rc.CWD:
            os.chdir('..')
