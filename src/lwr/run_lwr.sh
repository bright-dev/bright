#!/bin/sh
### PBS Settings
#PBS -N CHAR_lwr
#PBS -l ncpus=4,nodes=1:ppn=4,walltime=36:00:00
#PBS -k oe
#PBS -j oe
#PBS -W stagein=./lwr.i@nukestar01:/home/scopatz/runchar/lwr/lwr.i
#PBS -W stageout=./lwr.o@nukestar01:/home/scopatz/runchar/lwr/lwr.o
#PBS -W stageout=./lwr.s@nukestar01:/home/scopatz/runchar/lwr/lwr.s
#PBS -W stageout=./lwr.m@nukestar01:/home/scopatz/runchar/lwr/lwr.m
#PBS -W stageout=./lwr.r@nukestar01:/home/scopatz/runchar/lwr/lwr.r
#PBS -W stageout=./CHAR_lwr.o*@nukestar01:/home/scopatz/runchar/lwr/CHAR_lwr.o*
#PBS -M scopatz@gmail.com
#PBS -m abe

### Display the job context
echo ""
echo "Running on host" `hostname`
echo "Time is" `date`
echo "Directory is" `pwd`
echo "DATAPATH is ${DATAPATH}"
echo "The master node of this job is: $PBS_O_HOST"
NPROCS=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`
echo "This job is using $NPROCS CPU(s) on the following $NNODES node(s):"
echo "-----------------------"
uniq $PBS_NODEFILE | sort
echo "-----------------------"
echo ""

### Set MCNPX datapath variable
export DATAPATH=/usr/share/mcnpxv260/lib/

### Run MCNP with MPI
/usr/local/bin/mpiexec \
-machinefile $PBS_NODEFILE \
/usr/share/mcnpxv260/bin/mcnpx260 \
i=lwr.i \
o=lwr.o \
s=lwr.s \
m=lwr.m \
r=lwr.r   
