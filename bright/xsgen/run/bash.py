from __future__ import print_function

import os
import time
import subprocess


# Char Libraries
from runchar import RunChar
# from ..templates.run_script import template as run_script_template
from pyne.utils import message, failure


class Bash(RunChar):
    """A controller to run char in the Bourne-Again Shell."""

    def __init__(self, n_code, env):
        """Args:
            * n_code: a neutron transport model
            * env: the environment to execute the model in.
        """
        super(Bash, self).__init__(n_code, env)


    def make_run_script(self):
        """Makes a shell script that will execute char."""

        run_fill = {}

        run_fill["run_shell"] = "#!/bin/bash\n"
        run_fill["PBS_general_settings"] = ''
        run_fill['PBS_job_context']  = ''

        # Remote copy commands
        if self.env.LOCAL:
            run_fill['remote_put'] = ''
            run_fill['remote_get'] = ''
        else:
            run_fill['remote_put'] = ("scp -r {rc.user}@{rg}:{rc.dir} ~/tmpchar/\n"
                                      "cd ~/tmpchar/\n").format(rc=self.env.remote_connection,
                                                            rg=self.env.remote_gateway)
            run_fill['remote_get'] = ("cd ~\n"
                                      "mv ~/CHAR_*.o* ~/tmpchar/\n"
                                      "scp -r ~/tmpchar/* {rc.user}@{rg}:{rc.dir}\n"
                                      "rm -r ~/tmpchar/\n").format(rc=self.env.remote_connection,
                                                                   rg=self.env.remote_gateway)

        # Get transport specific values
        run_fill.update(self.n_code.run_script_fill_values())

        # Fill the template
        with open(self.env.run_script, 'w') as f:
            f.write(run_script_template.format(**run_fill))

        os.chmod(self.env.run_script, 0o755)


    #
    # Where to run methods
    #

    def run_locally(self):
        """Runs the transport calculation on the local machine."""
        t1 = time.time()
        subprocess.call("./{0}".format(self.env.run_script), shell=True)
        t2 = time.time()

        # Report times
        time_msg = "{0:.3G}".format((t2-t1)/60.0)
        self.env.logger.info("Transport executed in {0} minutes.".format(time_msg))
        if 0 < self.env.VERBOSE:
            print(message("\nTransport executed in {0} minutes.\n".format(time_msg)))


    def run_remotely(self):
        """Runs the transport calculation on a remote machine"""
        if 0 < self.env.VERBOSE:
            print(message("Copying files to remote server."))

        rs = self.env.run_script
        rc = self.env.remote_connection

        try:
            # Make remote directory, if it isn't already there
            rc.run("mkdir -p {rc.dir}".format(rc=rc))

            # Remove the current contents of the remote directory
            rc.run("rm -r {rc.dir}*".format(rc=rc))

            # Put all appropriate files in reomte dir
            rc.put(rs, rc.dir + rs)
            for inputfile in self.n_code.place_remote_files:
                print(inputfile)
                rc.put(inputfile, rc.dir + inputfile)

            # Run char
            rc.run("source /etc/profile; cd {rc.dir}; ./{rs} > run.log 2>&1 &".format(rc=rc,
                                                                                      rs=rs))
            if 0 < self.env.VERBOSE:
                print(message("Running transport code remotely."))

        except NameError:
            if 0 < self.env.VERBOSE:
                print(failure(remote_err_msg))

        # My work here is done
        raise SystemExit



    #
    # Process and file helper methods
    #

    def fetch(self):
        """Fetches files from remote server."""
        if 0 < self.env.VERBOSE:
            print(message("Fetching files from remote server."))

        try:
            for outputfile in self.n_code.fetch_remote_files:
                self.env.remote_connection.get(self.env.remote_connection.dir + outputfile, ".")
        except NameError:
            if self.env.VERBOSE:
                print(failure(remote_err_msg))

        raise SystemExit


    def pid(self):
        """Finds the process id for a currently running char instance."""
        if self.env.LOCAL:
            rflag = ''
            spout = subprocess.check_output("ps ux | grep {0}".format(self.n_code.run_str), shell=True)
        else:
            rflag = 'Remote '
            try:
                spout = subprocess.check_output("ssh {rc.user}@{rc.url} 'ps ux | grep {0}'".format(
                                                self.n_code.run_str, rc=self.env.remote_connection),
                                                shell=True)
            except NameError:
                if 0 < self.env.VERBOSE:
                    print(failure(remote_err_msg))

        spout = spout.split('\n')[:-1]
        if len(spout) < 1:
            if 0 < self.env.VERBOSE:
                print(message("{0}Process Not Running.".format(rflag)))
            raise SystemExit

        self.pid = spout[0].split()[1]
        self.prt = spout[0].split()[9]

        if 0 < self.env.VERBOSE:
            print(message("{0}Process ID: {1}".format(rflag, self.pid)))
            print(message("{0}Process Runtime: {1} min.".format(rflag, self.prt)))

        raise SystemExit


    def kill(self):
        """Kills the currently running char process."""
        # Grab the pid
        try:
            self.pid()
        except SystemExit:
            pass

        # Kill the process
        if self.env.LOCAL:
            rflag = ''
            spout = subprocess.check_output("kill {0}".format(self.pid), shell=True)
        else:
            rflag = 'Remote '
            spout = '\n'
            self.env.remote_connection.run("kill {0}".format(self.pid))

        spout = spout.split('\n')[:-1]

        if 0 < self.env.VERBOSE:
            print(spout)
            print(message("{0}Process Killed.".format(rflag)))

        raise SystemExit

