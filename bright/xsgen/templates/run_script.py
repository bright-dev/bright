template = """\
{run_shell}
{PBS_general_settings}

### Display the job context
echo 
echo "Running on host" `hostname`
echo "Time is" `date`
echo "Directory is" `pwd`
{transport_job_context}
{PBS_job_context}

{remote_put}
{run_commands}
{remote_get}
"""
