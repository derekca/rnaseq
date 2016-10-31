#!/bin/sh
# AUTHOR : Derek Caetano-Anolles
# EDITED : 2016.10
# USAGE  : sh ./06_sexing.sh <bam> <sexOUT>
# DESCR. : Reads bam files out of a 'bam' folder that should already exist. Finds
#          the Y-linked Eif2s3y gene, which determines sex. Exports a spreadsheet.

echo "./06_sexing.sh" # Echo to stdout.

#┌─────────────────────────────────────────────────────────────────────┐
#│ PASSED VARIABLES                                                    │
#└─────────────────────────────────────────────────────────────────────┘

bam=$1              # This bam output folder should already exist.
sexOUT=$2           # The output file for the sex determination.

#┌───────────┬──────────────────┬──────────────────────────────────────┐
#│ SAMTOOLS  │  0.1.19-44428cd  │  Determine sex.                      │
#└───────────┴──────────────────┴──────────────────────────────────────┘
# Finds the Y-linked Eif2s3y gene, which determines sex. Exports into "sexOUT".

for bamfiles in `ls $bam/*.bam`
do

    chrY="chrY:1010543-1028847"                   # The location of the gene being searched.
    bamNAMES=${bam##*/}                           # Names of the bamfiles.
    echo -ne ${bamNAMES%_*}"\t" > $sexOUT         # Output the file name to TSV file.
    samtools index ${bamfiles}                    # Create dictionary to search through.
    samtools view -c ${bamfiles} $chrY > $sexOUT  # Search index and count the number of
                                                  # times "chrY..." shows up. Write in TSV.
done


# LOGFILE
sh ./AL_logs.sh "log" "SamTools sexing" "$sexOUT"
