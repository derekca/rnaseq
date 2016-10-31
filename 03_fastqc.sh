#!/bin/sh
# AUTHOR : Derek Caetano-Anolles
# EDITED : 2016.10
# USAGE  : sh ./03_fastqc.sh <out> <QCout> <QCnames>
# DESCR. : Checks quality of RAW or TRIMMED fastq files.

echo "./03_fastqc.sh" # Echo to stdout.

#┌─────────────────────────────────────────────────────────────────────┐
#│ PASSED VARIABLES                                                    │
#└─────────────────────────────────────────────────────────────────────┘

out=$1                 # This is where the Fastq.gz (RAW or TRIM) files are located.
QCout=$2               # This is where the FastQC output goes.
QCnames=$3             # The naming convention for the files being QC'd.
                       # ie. "[1-9]*.fastq.gz"


#┌────────┬───────────┬────────────────────────────────────────────────┐
#│ FASTQC │  v0.11.4  │  Check RAW fastq quality.                      │
#└────────┴───────────┴────────────────────────────────────────────────┘

mkdir -p $QCout                  # Make directory

# Output FastQC results into [$QCout] dir, and run it on the RAW files,
# which are ID'd as following the naming pattern described in $rawNAMES.

fastqc -t 16 -o $QCout $out/$QCnames

# LOGFILE
sh ./AL_logs.sh "log" "FastQC_`basename $out`" "$QCout"