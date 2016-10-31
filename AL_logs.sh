#!/bin/sh
# AUTHOR : Derek Caetano-Anolles
# EDITED : 2016.10
# USAGE  : sh ./AL_logs.sh <module> <input1> <input2> <input2>
# DESCR. : This script is in charge of writing the logfile.
#          The logfile is written in a modular way, with different
#          logfile outputs being written based on which module is called.
#          For example, calling "ini" writes the job initiation output
#          at the start of the logfile. Here is an example of a full log:

#          BATCH JOB ID: 81841.darwin
#          Fri Oct 21 18:15:20 CEST 2016 | JOB INITIATED.
#          ──────────────────────────────────────────────
#          Fri Oct 21 18:15:20 CEST 2016 | PIPELINE cfw10_160530 INITIATED
#          Fri Oct 21 18:33:04 CEST 2016 |  [ ]SUCCESS bcl2fastq2.
#          Fri Oct 21 18:59:27 CEST 2016 |  [ ]SUCCESS FastQC_raw.
#          Fri Oct 21 22:15:56 CEST 2016 |  [ ]SUCCESS Trimmomatic.
#          Fri Oct 21 22:20:21 CEST 2016 |  [ ]SUCCESS FastQC_trim.
#          Sat Oct 22 03:12:53 CEST 2016 |  [ ]SUCCESS HISAT.
#          Sat Oct 22 03:12:53 CEST 2016 |  [ ]SUCCESS StringTie_fpkm.
#          Sat Oct 22 03:12:53 CEST 2016 |  [ ]SUCCESS StringTie_ctab.
#          Sat Oct 22 03:13:02 CEST 2016 |  [ ]SUCCESS SamTools sexing.
#          Sat Oct 22 03:13:02 CEST 2016 |   > 8h:57m:42s
#          ──────────────────────────────────────────────
#          Sat Oct 22 03:13:02 CEST 2016 | JOB COMPLETED.
#          EXIT CODE IS: 0

#┌─────────────────────────────────────────────────────────────────────┐
#│ PASSED VARIABLES                                                    │
#└─────────────────────────────────────────────────────────────────────┘

module=$1   # Describes a module below, to be run by the script.
input1=$2   # Name of the job that will be described in the log output.
input2=$3   # Generic variable that can be used by a module.
input3=$4   # Generic variable that can be used by a module (unused).


#┌─────────────────────────────────────────────────────────────────────┐
#│ LOG MODULES                                                         │
#└─────────────────────────────────────────────────────────────────────┘

#┌─────┐ MODULE: Make Logfile
#│     │ DESCR : << Happens every time this script is run. >>
#└─────┘ INPUT : Not applicable.

mkdir -p $HOME/std_logs                   # Make log folder (if it's not there).
logfile=$HOME/std_logs/"${PBS_JOBID}.log" # Save log under this name.
export logfile                            # Export log so all scripts can write in it.


#┌─────┐ MODULE: Job Initiation
#│ ini │ DESCR : Log a batch job inititiation, and start Log_Timer.
#└─────┘ INPUT : <$module> is 'ini'
#                <$input1> is NULL
#                <$input2> is NULL

if [ $module == "ini" ]; then

    echo "BATCH JOB ID: ${PBS_JOBID}" >> $logfile
    echo "`date` | JOB INITIATED." >> $logfile
    echo "──────────────────────────────────────────────" >> $logfile

fi


#┌─────┐ MODULE: Pipeline-Loop
#│ pll │ DESCR : Log the start of a pipeline loop.
#└─────┘ INPUT : <$module> is 'pll'
#                <$input1> is the output folder name (as an identifier).
#                <$input2> is NULL

if [ $module == "pll" ]; then

    echo "`date` | PIPELINE $input1 INITIATED" >> $logfile

fi


#┌─────┐ MODULE: Log Success
#│ log │ DESCR : Log a script's success, or clean up a failure.
#└─────┘ INPUT : <$module> is 'log'
#                <$input1> should be the software package that generates files.
#                <$input2> should be the name of the tool's output folder.

# If the folder has files in it, then log a "success" for $input1.
# If the folder is empty (which happens when there is an error/bug,
# then log a "failure" for $input1, and remove the folder.

if [ $module == "log" ]; then

    if [ "$(ls -A $input2)" ]; then

        echo "`date` |  [ ]SUCCESS $input1." >> $logfile
        else
        echo "`date` |  [*]FAILURE $input1." >> $logfile
        rm -r $input2

    fi

fi


#┌─────┐ MODULE: Pipeline-End
#│ ple │ DESCR : Log the end of a pipeline loop.
#└─────┘ INPUT : <$module> is 'ple'

if [ $module == "ple" ]; then

    echo "`date` | ALL PIPELINES ENDED." >> $logfile

fi


#┌─────┐ MODULE: Pipeline-Timer
#│ hrs │ DESCR : Log how much time the loop has taken.
#└─────┘ INPUT : <$module> is 'hrs'
#                <$input1> should be the difference between start/end time in seconds.

if [ $module == "hrs" ]; then

    # END_SCRIPT_TIME
    hrs=$(($input1/3600))          # Hrs in the timer.
    min=$(($input1%3600/60))       # Leftover min from the hrs.
    sec=$(($input1%60))            # Leftover sec from the min.

    echo "`date` |   > ${hrs}h:${min}m:${sec}s" >> $logfile

fi

#┌─────┐ MODULE: Job Completion
#│ end │ DESCR : Log a batch job completion.
#└─────┘ INPUT : <$module> is 'end'

if [ $module == "end" ]; then

    echo "──────────────────────────────────────────────" >> $logfile
    echo "`date` | JOB COMPLETED." >> $logfile
    echo "EXIT CODE IS: $?" >> $logfile

fi







