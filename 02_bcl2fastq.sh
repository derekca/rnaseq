#!/bin/sh
# AUTHOR : Derek Caetano-Anolles
# EDITED : 2016.10
# USAGE  : sh ./02_bcl2fastq.sh <raw> <runpath>
# DESCR. : Runs bcl2fastq on a sample.

echo "./02_bcl2fastq.sh" # Echo to stdout.


#┌─────────────────────────────────────────────────────────────────────┐
#│ PASSED VARIABLES                                                    │
#└─────────────────────────────────────────────────────────────────────┘

raw=$1                 # This is where the bcl2fastq2 files go.
runpath=$2             # The directory where all SEQ DATA is located.

#┌────────────┬──────────────┬─────────────────────────────────────────┐
#│ BCL2FASTQ  │  v2.17.1.14  │  Converts bcl to fastq.                 │
#└────────────┴──────────────┴─────────────────────────────────────────┘

mkdir -p $raw                       # Make directory

bcl2fastq -R "$runpath" -o "$raw"   # Run bcl2fastq
                                    # [-R] with NextSeq runpath, and
                                    # [-o] set output path to $raw

# LOGFILE
sh ./AL_logs.sh "log" "bcl2fastq2" "$raw"










