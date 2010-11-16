import os

######################
### CHAR Libraries ###
######################
from char import *


def make_run_script(n_transporter, runflag, localflag=True):
    run_fill = {}

    if runflag in ["PBS"]:
        run_fill["Run_Shell"] = "#!/bin/sh"

        run_fill["PBS_General_Settings"] = "### PBS Settings\n"
        run_fill["PBS_General_Settings"] = "#PBS -N CHAR_{0}\n".format(reactor)

        run_fill["PBS_General_Settings"] = "#PBS -l ncpus={0}".format(NumberCPUs)
        run_fill["PBS_General_Settings"] = ",nodes={0}".format(NumberCPUs/CPUsPerNode)
        run_fill["PBS_General_Settings"] = ":ppn={0}".format(CPUsPerNode)

        run_fill["PBS_General_Settings"] = "#PBS -k oe\n"
        run_fill["PBS_General_Settings"] = "#PBS -j oe\n"

        run_fill["PBS_General_Settings"] = "#PBS -M {0}\n".format(email)
        run_fill["PBS_General_Settings"] = "#PBS -m abe\n".format(email)
    else:
        run_fill["Run_Shell"] = "#!/bin/bash\n"
        run_fill["PBS_General_Settings"] = ''

    # Set PBS_Job_Context
    if runflag in ["PBS"]:
        run_fill['PBS_Job_Context']  = "echo \"The master node of this job is: $PBS_O_HOST\"\n"
        run_fill['PBS_Job_Context'] += "NPROCS=`wc -l < $PBS_NODEFILE`\n"
        run_fill['PBS_Job_Context'] += "NNODES=`uniq $PBS_NODEFILE | wc -l`\n"
        run_fill['PBS_Job_Context'] += "echo \"This job is using $NPROCS CPU(s) on the following $NNODES node(s):\"\n"
        run_fill['PBS_Job_Context'] += "echo \"-----------------------\"\n"
        run_fill['PBS_Job_Context'] += "uniq $PBS_NODEFILE | sort\n"
        run_fill['PBS_Job_Context'] += "echo \"-----------------------\"\n"
        run_fill['PBS_Job_Context'] += "echo \"\"\n"
    else:
        run_fill['PBS_Job_Context']  = ''

    run_fill.update(n_transporter.run_script_fill_values(runflag))

    # Fill the template
    with open('../templates/run_script.sh.template', 'r') as f:
        run_script_template = f.read()

    with open(runscript, 'w') as f:
        f.write(run_script_template.format(**run_fill))

    os.chmod(runscript, 0755)

    return

def run_transport_local(runflag):
    """Runs the transport calculation on the local machine."""
    t1 = time.time()
    if runflag == "PBS":
        subprocess.call("qsub {0}".format(runscript), shell=True)
    else:
        subprocess.call("./{0}".format(runscript), shell=True)
    t2 = time.time()
    if 0 < verbosity:
        print()
        print(message("Transport executed in {0:time} minutes.", "{0:.3G}".format((t2-t1)/60.0) ))
        print()
    return

def run_transport_remote(runflag):
    """Runs the transport calculation on a remote machine"""
    try:
        if 0 < verbosity:
            print(message("Copying files to remote server."))

        RemoteConnection.run("mkdir -p {rc.RemoteDir}".format(rc=RemoteConnection)) 
        RemoteConnection.run("rm {rc.RemoteDir}*".format(rc=RemoteConnection))
        RemoteConnection.put(runscript,  RemoteConnection.RemoteDir + runscript)
        for inputfile in n_transporter.place_remote_files:
            RemoteConnection.put(inputfile,  RemoteConnection.RemoteDir + inputfile)

        if runflag == "PBS":
            RemoteConnection.run("source /etc/profile; cd {rc.RemoteDir}; qsub {rs} > runlog.txt 2>&1 &".format(rc=RemoteConnection, rs=runscript))
        else:
            RemoteConnection.run("source /etc/profile; cd {rc.RemoteDir}; ./{rs} > runlog.txt 2>&1 &".format(rc=RemoteConnection, rs=runscript))

        if 0 < verbosity:
            print(message("Running transport code remotely."))
        raise SystemExit
    except NameError:
        if 0 < verbosity:
            print(failure("Host, username, password, or directory not properly specified for remote machine."))
            print(failure("Please edit defchar to include RemoteURL, RemoteUser, RemotePass, and RemoteDir."))
            print(failure("Nothing to do, quiting."))
        raise SystemExit
    return

def fetch_remote_files():
    """Fetches files from remote server."""
    try:
        if 0 < verbosity:
            print(message("Fetching files from remote server."))
        for outputfile in n_transporter.fetch_remote_files:
            metasci.SafeRemove(outputfile)
        for outputfile in n_transporter.fetch_remote_files:
            RemoteConnection.get(RemoteConnection.RemoteDir + outputfile, ".")
    except NameError:
        if 0 < verbosity:
            print(message("Host, username, password, or directory not properly specified for remote machine."))
            print(message("Please edit defchar to include RemoteURL, RemoteUser, RemotePass, and RemoteDir."))
            print(message("Nothing to do, quiting."))
        raise SystemExit
    return

def find_pid_local(BoolKill = False):
    """Finds (and kills?) the Local Transport Run Process"""
    sp = subprocess.Popen("ps ux | grep {0}".format(n_transporter.run_str), stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, shell=True) 
    spout, sperr = sp.communicate() 
    spout = spout.split('\n')[:-1]
    pid =  spout[0].split()[1]
    prt =  spout[0].split()[9]
    if 0 < verbosity:
        print(message("Process ID: {0}".format(pid)))
        print(message("Process Runtime: {0:time} min.".format(prt)))
    del sp, spout, sperr

    if BoolKill:
        sp = subprocess.Popen("kill {0}".format(pid), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True) 
        spout, sperr = sp.communicate() 
        spout = spout.split('\n')[:-1]
        if 0 < verbosity:
            print(spout)
            print(message("Process Killed."))
        del sp, spout, sperr
        raise SystemExit
    return

def find_pid_remote(runflag, BoolKill = False):
    """Finds (and kills?) the Remote Transport Run Process"""
    try:
        if runflag in ["PBS"]:
            rsp = subprocess.Popen("ssh {rc.RemoteUser}@{rc.RemoteURL} \"qstat -u {rc.RemoteUser}\"".format(rc=RemoteConnection), 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True) 
        else:
            rsp = subprocess.Popen("ssh {rc.RemoteUser}@{rc.RemoteURL} \"ps ux | grep {0}\"".format(n_transporter.run_str, 
                rc=RemoteConnection), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True) 
    except NameError:
        if 0 < verbosity:
            print(failure("Host, username, password, or directory not properly specified for remote machine."))
            print(failure("Please edit defchar to include RemoteURL, RemoteUser, RemotePass, and RemoteDir."))
            print(failure("Nothing to do, quiting."))
        raise SystemExit

    #grab and parse the remote process info
    spout, sperr = rsp.communicate() 
    if runflag in ["PBS"]:
        spl = spout.split('\n')
        if len(spl) < 2:
            if 0 < verbosity:
                print(message("Remote Process Not Running."))
            raise SystemExit
            return
        spll = spl[-2].split()
        pid = spll[0].partition('.')[0]
        prt = spll[-1]
    else:
        spl = spout.split('\n')
        if len(spl) < 1:
            if 0 < verbosity:
                print(message("Remote Process Not Running."))
            raise SystemExit
            return
        spll = spl[:-1]
        pid =  spll[0].split()[1]
        prt =  spll[0].split()[9]

    #Print the process info
    if 0 < verbosity:
        print(message("Remote Process ID: {0}".format(pid)))
        print(message("Remote Process Runtime: {0:time} min.".format(prt)))
        if runflag in ["PBS"]:
            print()
            print(spout)

    #Kill the remote process if required...
    if BoolKill:
        if runflag in ["PBS"]:
            RemoteConnection.run("qdel {0}".format(pid)) 
            RemoteConnection.run("cluster-kill mcnpx260") 
        else:
            RemoteConnection.run("kill {0}".format(pid)) 

        if 0 < verbosity:
            print(message("Remote Process Killed."))

    raise SystemExit
    return
