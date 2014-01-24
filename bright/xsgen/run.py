from bright.xsgen.plugins import Plugin

class XDressPlugin(Plugin):
    def update_argparser(self, parser):
        parser.add_argument("-r", "--run", action="store_true", dest="RUN_TRANSPORT",
            help="Run the transport calculation.")

        parser.add_argument("-k", "--kill", action="store_true", dest="KILL_TRANSPORT",
            help="Kills the current transport run. Sets -p.")

        parser.add_argument("-p", "--pid", action="store_true", dest="PID",
            help="Finds the process identification number of a current transport run. Sets -d.")

        parser.add_argument("-f", "--fetch", action="store_true", dest="FETCH_FILES",
            help="Fetches files from the remote server. Does not run transport, even if -r is set. Automatically sets -s.")

        parser.add_argument("-b", "--burnup",  action="store_true", dest="RUN_BURNUP",
            help="Run the burnup calculation.")

        parser.add_argument("-x", "--xs",  action="store_true", dest="RUN_XS_GEN", 
            help="Run the cross-section generation calculation.")

        parser.add_argument("-m", "--delta-mass",  action="store_true", dest="RUN_DELTAM", 
            help="Run the initial nuclide sensitivity calculation.")

        parser.add_argument("-a", "--analyze", action="store_true", dest="RUN_ANALYSIS", 
            help="Run analysis on database.")

    def setup(self, rc):
        if rc.KILL_TRANSPORT:
            rc.PID = True                      #Ensures that the PID is found in order that it mak be killed.

        if rc.PID:
            rc.RUN_XS_GEN = False
            rc.RUN_DELTAM = False
            rc.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.

        if rc.FETCH_FILES:
            rc.RUN_XS_GEN = False
            rc.RUN_DELTAM = False
            rc.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.
            rc.LOCAL = False                   #Ensures that ssh package is loaded.

    def execute(self, rc):
        if rc.RUN_TRANSPORT:
            # Run Transport code
            rc.runchar.make_run_script()

            if rc.LOCAL:
                rc.runchar.run_locally()
            else:
                rc.runchar.run_remotely()

        elif rc.RUN_ANALYSIS:
            n_code.analyze_deltam()

        elif rc.RUN_BURNUP or rc.RUN_XS_GEN or rc.RUN_DELTAM:
            # Make tranumatrion libraries by executing the as a separate step from 
            # the cross-section generation
            if rc.RUN_BURNUP:
                rc.runchar.burnup(idx)

            # Make Cross-sections as a separate step from the burnup calculation
            if rc.RUN_XS_GEN:
                rc.runchar.xs_gen(idx, nucs)

            # Run initial nuclide sensitivity calculation
            if rc.RUN_DELTAM:
                n_code.run_deltam_pert(idx, ihm_isos, sidx)

        elif rc.FETCH_FILES:
            # Fetches files from remote server
            rc.runchar.fetch()
        elif rc.PID:
            # Finds and prints the PID of CHAR
            rc.runchar.pid()
        elif rc.KILL_TRANSPORT:
            # Finds and kills CHAR
            rc.runchar.kill()
