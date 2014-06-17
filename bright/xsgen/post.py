from bright.xsgen.plugins import Plugin

class XSGenPlugin(Plugin):
    def execute(self, rc):
        # Clean up
        if not rc.CWD:
            os.chdir('..')
