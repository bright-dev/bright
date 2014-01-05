from plugins import Plugin

class XDressPlugin(Plugin):

    defaultrc = {"run_switch": {'': Bash,
                                'BASH': Bash,
                                'bash': Bash,
                                'PBS': Pbs,
                                'pbs': Pbs,
                                'Torque': Pbs,
                                'torque': Pbs}}
    def update_argparser(self, parser):

        parser.add_option("-i", "--input", action="store_true", dest="MAKE_INPUT",
                          help="Makes the transport calculation input deck.")

        # default=False needs to go in default rc
        parser.add_option("-r", "--run", action="store_true", dest="RUN_TRANSPORT",
            help="Run the transport calculation.")

        parser.add_option("-d", "--dry-run", action="store_false", dest="RUN_TRANSPORT",
            help="Dry Run. Do NOT run the transport calculation.")

        parser.add_option("-k", "--kill", action="store_true", dest="KILL_TRANSPORT", 
            default=False, help="Kills the current transport run. Sets -p.")


    def execute(self, rc=self.defaultrc):
        runner = rc["run_switch"]["scheduler"]
        runchar = runner(n_code, rc)
        if (options.MAKE_INPUT) and (not options.FETCH_FILES) and (not options.PID):
            runchar.init_h5()

