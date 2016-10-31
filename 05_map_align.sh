#!/bin/sh
# AUTHOR : Derek Caetano-Anolles
# EDITED : 2016.10
# USAGE  : sh ./05_map_align.sh <trim> <bam> <fpkm> <ctab> <hisatidx> <refannot>
# DESCR. : Uses HISAT and StringTie to map/align reads that are found in a
#          "trim" folder.

echo "./05_map_align.sh" # Echo to stdout.

#┌─────────────────────────────────────────────────────────────────────┐
#│ PASSED VARIABLES                                                    │
#└─────────────────────────────────────────────────────────────────────┘

trim=$1      # This is where the TRIMMED reads should already exist.
bam=$2       # This is where the BAM files go.
fpkm=$3      # This is where the FPKM files go.
ctab=$4      # This is where the CTAB files go.

hisatidx=$5  # The directory for the reference genome file.
refannot=$6  # The directory for the reference gene annotation file.

#┌─────────────────────────────────────────────────────────────────────┐
#│ MAKE DIRECTORIES                                                    │
#└─────────────────────────────────────────────────────────────────────┘

mkdir -p $bam
mkdir -p $fpkm
mkdir -p $ctab

#┌─────────────────────────────────────────────────────────────────────┐
#│ MAP AND ALIGN READS                                                 │
#└─────────────────────────────────────────────────────────────────────┘

# For all the trimmed "L001_R1" files, run HISAT and StringTie.
for fastq1 in `ls $trim/*L001_R1*`
do

    #┌───────────┬──────────────────┬──────────────────────────────────┐
    #│ HISAT     │  v2.0.1-beta     │  Map reads.                      │
    #└───────────┴──────────────────┴──────────────────────────────────┘
    #┌───────────┬──────────────────┬──────────────────────────────────┐
    #│ SAMtools  │  0.1.19-44428cd  │  Sort alignments.                │
    #└───────────┴──────────────────┴──────────────────────────────────┘

    # Naming convention for indexing.
    fastq2=${fastq1/R1/R2}  # fastq1 is all R1 files, fastq2 is all R2 files.
    tmp1=${fastq1##*/}      # tmp1 is all the trim files with L001_R1.
    tmp2=${tmp1#*_}         # tmp2 is all the L001_R1
    output=${tmp2%_L001*}   # output is the name of the outputed file (also used for StringTie).

    # Run HISAT2.
    hisat2 -x $hisatidx \
        -1 ${fastq1},${fastq1/L001/L002},${fastq1/L001/L003},${fastq1/L001/L004} \
        -2 ${fastq2},${fastq2/L001/L002},${fastq2/L001/L003},${fastq2/L001/L004} \
        --rna-strandness RF \
        --dta \
        -p 16 | \
        samtools view -Su - | \
        samtools sort - $bam/${output}_sorted

    #┌────────────┬──────────┬─────────────────────────────────────────┐
    #│ STRINGTIE  │  v1.2.1  │  Sequence assembly, estimate abundances │
    #└────────────┴──────────┴─────────────────────────────────────────┘
    # Run StringTie using the bam files located at the given directory.
    # Set names of GTF output [-o]

    stringtie $bam/${output}_sorted.bam \
        -o "$fpkm"/${output}.gtf \
        -p 16 \
        -G "$refannot" \
        -b "$ctab"/${output} \
        -e

done

#!  # Remove trim directory to conserve space.
    # rm -r $trim


# LOGFILE
sh ./AL_logs.sh "log" "HISAT" "$bam"
sh ./AL_logs.sh "log" "StringTie_fpkm" "$fpkm"
sh ./AL_logs.sh "log" "StringTie_ctab" "$ctab"




