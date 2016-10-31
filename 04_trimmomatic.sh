#!/bin/sh
# AUTHOR : Derek Caetano-Anolles
# EDITED : 2016.10
# USAGE  : sh ./04_trimmomatic.sh <raw> <trim> <adapters> <trimNAMES>
# DESCR. : Runs Trimmomatic on files in the "raw" folder.

echo "./04_trimmomatic.sh" # Echo to stdout.

#┌─────────────────────────────────────────────────────────────────────┐
#│ PASSED VARIABLES                                                    │
#└─────────────────────────────────────────────────────────────────────┘

raw=$1              # This is where the RAW reads should already exist.
trim=$2             # This is where the TRIMMED reads go.

adapters=$3         # This is the filepath for the ADAPTERS FILE.
trimNAMES=$4        # This is the naming convention for the RAW files
                    # being trimmed, ie. "[1-9]*.fastq.gz"

#┌─────────────────────────────────────────────────────────────────────┐
#│ MAKE DIRECTORIES                                                    │
#└─────────────────────────────────────────────────────────────────────┘

mkdir -p $trim

#┌──────────────┬─────────┬────────────────────────────────────────────┐
#│ TRIMMOMATIC  │  v0.35  │  Trim reads.                               │
#└──────────────┴─────────┴────────────────────────────────────────────┘

# Directory of Trimmomatic JAR
DIRtrimmomatic="$HOME/../../file/path/to/trimmomatic-0.35.jar"

# Run Trimmomatic JAR on the raw Fastq files, into the trimmed folder.
for fastq1 in `ls $raw/$trimNAMES`
do



    # Define variables for basein <inputR1><inputR2>
    # fastq1 is just the R1 values, fastq2 is the R2 values.
    fastq2=${fastq1/R1/R2}
    # Define prefix for the baseout name
    tmp=${fastq1##*/}
    output=${tmp%_R1*}

    # Run Trimmomatic.
    java -jar "$DIRtrimmomatic" \
        PE \
        -threads 16 \
        -phred33 \
        ${fastq1} \
        ${fastq2} \
        $trim/${output}_R1p.fastq.gz \
        $trim/${output}_R1u.fastq.gz \
        $trim/${output}_R2p.fastq.gz \
        $trim/${output}_R2u.fastq.gz \
        ILLUMINACLIP:"$adapters":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36

    # Remove the unpaired "R1u" and "R2u" files.
    rm $trim/${output}_R1u.fastq.gz \
       $trim/${output}_R2u.fastq.gz



done

# LOGFILE
sh ./AL_logs.sh "log" "Trimmomatic" "$trim"



