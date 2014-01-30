from bright.xsgen.plugins import Plugin

class XDressPlugin(Plugin):
    def update_argparser(self, parser):
        parser.add_argument("-i", "--input", action="store_true", dest="MAKE_INPUT",
            help="Makes the transport calculation input deck.")

    def execute(self, rc):
        # Get the run controller
        runner = run_switch[rc.scheduler]
        runchar = runner(n_code, rc)

        run_switch = {'': Bash,
                      'BASH': Bash,
                      'bash': Bash,
                      'PBS': Pbs,
                      'pbs': Pbs,
                      'Torque': Pbs,
                      'torque': Pbs,
                      }

        if (rc.MAKE_INPUT) and (not rc.FETCH_FILES) and (not rc.PID):
            runchar.init_h5()
