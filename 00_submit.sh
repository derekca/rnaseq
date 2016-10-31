#!/bin/sh
# AUTHOR : Derek Caetano-Anolles
# EDITED : 2016.10
# USAGE  : qsub ./00_submit.sh <suffix>
# DESCR. : Generic script for running any batch job to a queue. Can be
#          modified with script names and passed arguments in the "RUN
#          SCRIPT" section.

#┌─────────────────────────────────────────────────────────────────────┐
#│ PORTABLE BATCH SYSTEM (#PBS) INSTRUCTIONS                           │
#└─────────────────────────────────────────────────────────────────────┘

#PBS -S /bin/sh
#PBS -q serial
#PBS -V
#PBS -l mem=16GB
#PBS -l walltime=144:00:00
#PBS -l nodes=1:ppn=1
#PBS -o stdout.log
#PBS -e stderr.log
#PBS -N cfw_pipe

#┌─────────────────────────────────────────────────────────────────────┐
#│ PASSED VARIABLES                                                    │
#└─────────────────────────────────────────────────────────────────────┘

suffix="$suffix"     # If there is a job-suffix then it will get passed along.
                     # It should come from the job-splitting script that
                     # qsubs the split-jobs, and it needs to be global.

#┌─────────────────────────────────────────────────────────────────────┐
#│ GENERATE LOG DATA                                                   │
#└─────────────────────────────────────────────────────────────────────┘

cd $PBS_O_WORKDIR                  # Enter $PBS working directory.
rm stdout.log                      # Remove previous stdout.log file.
rm stderr.log                      # Remove previous stderr.log file.

# START_STDOUT - log environmental data.
echo "PBS JOBID : $PBS_JOBID"
echo "JOBNAME   : $PBS_JOBNAME"
echo "JOBSUFFIX : $suffix"
echo "DATE      : `date`"
echo "WORKDIR   : $PBS_O_WORKDIR"
echo "PWD IS    : `pwd`"
echo "HOSTNAME  : `hostname`"
echo "PPN USED  :" `cat $PBS_NODEFILE`
echo "INITIATE  : `date`"
echo "┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄"


# START_SCRIPT_TIMER
time0=$SECONDS  # Save/export the start time of the script, so it can
export time0    # be used for calculating how long the script runs.

# START_LOGFILE
sh ./AL_logs.sh "ini"

# SPECIFY STOUT REDIRECT
# Add '>$logerr 2>&1' after running a script in order to
# redirect all of that script's outputs into this file.
logerr="$HOME/std_logs/${PBS_JOBID}.err.log"



#╔═════════════════╗┌──────────────────────────────────────────────────┐
#║  RUN SCRIPT(S)  ║│ Run the script filepaths in this space.          │
#╚═════════════════╝└──────────────────────────────────────────────────┘
#├─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─┤





sh ./01_PIPELINE.sh "$suffix" >$logerr 2>&1






#├─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─┤

#┌─────────────────────────────────────────────────────────────────────┐
#│ END SCRIPT                                                          │
#└─────────────────────────────────────────────────────────────────────┘

# END_SCRIPT_TIME
timeX=$(($SECONDS - $time0)) # How long this script took to run in sec.
hrs=$(($timeX/3600))         # Hrs in the timer.
min=$(($timeX%3600/60))      # Leftover min from the hrs.
sec=$(($timeX%60))           # Leftover sec from the min.

# END_LOGFILE
sh ./AL_logs.sh "end"

# END_STDOUT - Add confirmation of completion.
printf "┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄ \n"
printf "COMPLETED : `date` \n"
printf 'TOTALTIME : %02d:%02d:%02d \n' $hrs $min $sec
printf "EXIT CODE : $? \n"
if [ $? = 0 ]; then
echo '         '
echo '  ()-()  ' # Mouse appears
echo '   \"/   ' # when the exit
echo '    `    ' # code is 0. :)
fi

# Wait for background jobs to complete.
wait








