template = """\
{Run_Shell}
{PBS_General_Settings}{PBS_Walltime}{PBS_Stagein_Settings}{PBS_Stageout_Settings}

### Display the job context
echo 
echo "Running on host" `hostname`
echo "Time is" `date`
echo "Directory is" `pwd`
{Transport_Job_Context}
{PBS_Job_Context}

{Run_Commands}"""
