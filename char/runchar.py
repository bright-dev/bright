from __future__ import print_function

import os
import time
import subprocess
from math import ceil

from metasci.colortext import message, failure

######################
### CHAR Libraries ###
######################
from char import defchar
from templates.run_script import template as run_script_template

remote_err_msg = ("Host, username, password, or directory not properly specified for remote machine.\n"
                  "Please edit defchar to include remote_url, remote_user, and remote_dir.\n"
                  "Nothing to do, quiting.")


def make_run_script(n_transporter):
    global defchar
    from char import defchar

    run_fill = {}

    if defchar.scheduler in ["PBS"]:
        nodes = int(ceil( float(defchar.number_cpus) / defchar.cpus_per_node))

        run_fill["run_shell"] = "#!/bin/sh"

        run_fill["PBS_general_settings"]  = ""
        run_fill["PBS_general_settings"] += "### PBS Settings\n"
        run_fill["PBS_general_settings"] += "#PBS -N CHAR_{0}\n".format(defchar.reactor)

        run_fill["PBS_general_settings"] += "#PBS -l ncpus={0}".format(defchar.number_cpus)
        run_fill["PBS_general_settings"] += ",nodes={0}".format(nodes)
        run_fill["PBS_general_settings"] += ":ppn={0}".format(defchar.cpus_per_node) 
        run_fill["PBS_general_settings"] += ",walltime={0:02G}:00:00,pmem=3gb\n".format(
                                             n_transporter.run_script_walltime()) 

        run_fill["PBS_general_settings"] += "#PBS -k oe\n"
        run_fill["PBS_general_settings"] += "#PBS -j oe\n"

        run_fill["PBS_general_settings"] += "#PBS -M {0}\n".format(defchar.email)
        run_fill["PBS_general_settings"] += "#PBS -m abe\n".format(defchar.email)
    else:
        run_fill["run_shell"] = "#!/bin/bash\n"
        run_fill["PBS_general_settings"] = ''

    # Set PBS_Job_Context
    if defchar.scheduler in ["PBS"]:
        run_fill['PBS_job_context']  = "echo \"The master node of this job is: $PBS_O_HOST\"\n"
        run_fill['PBS_job_context'] += "NPROCS=`wc -l < $PBS_NODEFILE`\n"
        run_fill['PBS_job_context'] += "NNODES=`uniq $PBS_NODEFILE | wc -l`\n"
        run_fill['PBS_job_context'] += "echo \"This job is using $NPROCS CPU(s) on the following $NNODES node(s):\"\n"
        run_fill['PBS_job_context'] += "echo \"-----------------------\"\n"
        run_fill['PBS_job_context'] += "uniq $PBS_NODEFILE | sort\n"
        run_fill['PBS_job_context'] += "echo \"-----------------------\"\n"
        run_fill['PBS_job_context'] += "echo \"\"\n"
    else:
        run_fill['PBS_job_context']  = ''

    # Remote copy commands
    if defchar.options.LOCAL:
        run_fill['remote_put'] = ''
        run_fill['remote_get'] = ''
    else:
        run_fill['remote_put'] = ("scp -r {rc.user}@{rg}:{rc.dir} ~/tmpchar/\n"
                                  "cd ~/tmpchar/\n").format(rc=defchar.remote_connection, 
                                                            rg=defchar.remote_gateway)
        run_fill['remote_get'] = ("cd ~\n"
                                  "mv ~/CHAR_*.o* ~/tmpchar/\n"
                                  "scp -r ~/tmpchar/* {rc.user}@{rg}:{rc.dir}\n"
                                  "rm -r ~/tmpchar/\n").format(rc=defchar.remote_connection, 
                                                               rg=defchar.remote_gateway)

    # Get transport specific values
    run_fill.update(n_transporter.run_script_fill_values())

    # Fill the template
    with open(defchar.run_script, 'w') as f:
        f.write(run_script_template.format(**run_fill))

    os.chmod(defchar.run_script, 0755)

    return

def run_transport_local():
    """Runs the transport calculation on the local machine."""
    global defchar
    from char import defchar

    t1 = time.time()
    if defchar.scheduler == "PBS":
        subprocess.call("qsub {0}".format(defchar.run_script), shell=True)
    else:
        subprocess.call("./{0}".format(defchar.run_script), shell=True)
    t2 = time.time()

    # Report times
    time_msg = "{0:.3G}".format((t2-t1)/60.0)
    defchar.logger.info("Transport executed in {0} minutes.".format(time_msg))
    if 0 < defchar.verbosity:
        print()
        print(message("Transport executed in {0:time} minutes.", time_msg ))
        print()

    return

def run_transport_remote():
    """Runs the transport calculation on a remote machine"""
    global defchar
    from char import defchar

    try:
        if 0 < defchar.verbosity:
            print(message("Copying files to remote server."))

        # Make remote directory, if it isn't already there
        defchar.remote_connection.run("mkdir -p {rc.dir}".format(rc=defchar.remote_connection))

        # Remove the current contents of the remote directory        
        defchar.remote_connection.run("rm -r {rc.dir}*".format(rc=defchar.remote_connection))

        # Put all appropriate files in reomte dir
        defchar.remote_connection.put(defchar.run_script,  defchar.remote_connection.dir + defchar.run_script)
        for inputfile in defchar.n_transporter.place_remote_files:
            print(inputfile)
            defchar.remote_connection.put(inputfile,  defchar.remote_connection.dir + inputfile)

        if defchar.scheduler == "PBS":
            defchar.remote_connection.run("source /etc/profile; cd {rc.dir}; qsub {rs} > run.log 2>&1 &".format(
                                          rc=defchar.remote_connection, rs=defchar.run_script))
        else:
            defchar.remote_connection.run("source /etc/profile; cd {rc.dir}; ./{rs} > run.log 2>&1 &".format(
                                          rc=defchar.remote_connection, rs=defchar.run_script))

        if 0 < defchar.verbosity:
            print(message("Running transport code remotely."))

    except NameError:
        if 0 < defchar.verbosity:
            print(failure(remote_err_msg))

    # My work here is done
    raise SystemExit

def fetch_remote_files():
    """Fetches files from remote server."""
    global defchar
    from char import defchar

    if 0 < defchar.verbosity:
        print(message("Fetching files from remote server."))

    try:
        for outputfile in defchar.n_transporter.fetch_remote_files:
            defchar.remote_connection.get(defchar.remote_connection.dir + outputfile, ".")

    except NameError:
        if 0 < defchar.verbosity:
            print(failure(remote_err_msg))

        raise SystemExit

def find_pid_local():
    """Finds (and kills?) the Local Transport Run Process"""
    global defchar
    from char import defchar

    sp = subprocess.Popen("ps ux | grep {0}".format(n_transporter.run_str), stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE, shell=True) 
    spout, sperr = sp.communicate() 
    spout = spout.split('\n')[:-1]
    pid =  spout[0].split()[1]
    prt =  spout[0].split()[9]

    if 0 < defchar.verbosity:
        print(message("Process ID: {0}".format(pid)))
        print(message("Process Runtime: {0:time} min.".format(prt)))

    del sp, spout, sperr

    if defchar.options.KILL_TRANSPORT:
        sp = subprocess.Popen("kill {0}".format(pid), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        spout, sperr = sp.communicate() 
        spout = spout.split('\n')[:-1]

        if 0 < defchar.verbosity:
            print(spout)
            print(message("Process Killed."))

        del sp, spout, sperr
        raise SystemExit

def find_pid_remote():
    """Finds (and kills?) the Remote Transport Run Process"""
    global defchar
    from char import defchar

    try:
        if defchar.scheduler in ["PBS"]:
            rsp = subprocess.Popen("ssh {rc.user}@{rc.url} 'qstat -u {rc.user}'".format(
                                   rc=defchar.remote_connection), stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, shell=True) 
        else:
            rsp = subprocess.Popen("ssh {rc.user}@{rc.url} 'ps ux | grep {0}'".format(n_transporter.run_str, 
                                   rc=defchar.remote_connection), stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, shell=True) 
    except NameError:
        if 0 < defchar.verbosity:
            print(failure(remote_err_msg))

        raise SystemExit

    #grab and parse the remote process info
    spout, sperr = rsp.communicate() 
    if defchar.scheduler in ["PBS"]:
        spl = spout.split('\n')

        if len(spl) < 2:
            if 0 < defchar.verbosity:
                print(message("Remote Process Not Running."))

            raise SystemExit

        spll = spl[-2].split()
        pid = spll[0].partition('.')[0]
        prt = spll[-1]

    else:
        spl = spout.split('\n')
        if len(spl) < 1:

            if 0 < defchar.verbosity:
                print(message("Remote Process Not Running."))

            raise SystemExit

        spll = spl[:-1]
        pid =  spll[0].split()[1]
        prt =  spll[0].split()[9]

    #Print the process info
    if 0 < defchar.verbosity:
        print(message("Remote Process ID: {0}".format(pid)))
        print(message("Remote Process Runtime: {0:time} min.", prt))

        if defchar.scheduler in ["PBS"]:
            print(spout)

    #Kill the remote process if required...
    if defchar.options.KILL_TRANSPORT:
        if defchar.scheduler in ["PBS"]:
            defchar.remote_connection.run("qdel {0}".format(pid)) 

            if ('serpent' in defchar.transport_code) and ('mcnp' not in defchar.transport_code):
                defchar.remote_connection.run("cluster-kill sss-dev")
            elif ('mcnp' in defchar.transport_code) and ('serpent' not in defchar.transport_code):
                defchar.remote_connection.run("cluster-kill mcnpx260")

        else:
            defchar.remote_connection.run("kill {0}".format(pid)) 

        if 0 < defchar.verbosity:
            print(message("Remote Process Killed."))

    raise SystemExit
