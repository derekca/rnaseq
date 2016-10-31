#!/bin/sh
# AUTHOR : Derek Caetano-Anolles
# EDITED : 2016.10
# USAGE  : sh ./01_pipeline.sh <suffix>
# DESCR. : Runs script for project data. Specify the inputs for all
#          subscripts and it will carry over to the others.

echo "./01_pipeline.sh" # Echo to stdout.

#┌─────────────────────────────────────────────────────────────────────┐
#│ PASSED VARIABLES                                                    │
#└─────────────────────────────────────────────────────────────────────┘

suffix=$1       # If there is a suffix then the script will try and find an input file
                # named 'input${suffix}.sh' and run the inputs inside it.
                # Otherwise, it will just run 'suffix.sh'. This means that you can use
                # 00_submit.sh to run multiple samples through the pipeline at once, or
                # you can use 00_qsub_split the jobs in "index.sh" first, before running
                # 00_submit.sh on each split job individually.

#┌─────────────────────────────────────────────────────────────────────┐
#│ SCRIPT PIPELINE                                                     │
#└─────────────────────────────────────────────────────────────────────┘
# Read each line of $inputfile and perform the following scripts on the
# samples described in each line in that document, line-by-line.

inputfile="$HOME/import/input${suffix}.sh"      # CSV with two columns - $C1exportdir, $C2sampledir.
                                                # C1exportdir is the export folder name.
                                                # C2sampledir is the folder name where the samples are.

while read C1exportdir C2sampledir              # Read the 2 columns from the $inputfile in this loop.
do
    # LOGFILE
    loop0=$SECONDS                              # Loop timer starts here.
    sh ./AL_logs.sh "pll" "$C1exportdir"        # Log the start of the current loop into the logfile.

    #┌─────────────────┐
    #│  CONFIGURATION  ├─ - - - - - - - - - - - - - - - - - - - - - - -
    #└─────────────────┘

    # READ FROM INPUT.SH, DON'T EDIT
    C1exportdir="$C1exportdir"                  # Comes from the inputfile being read right now.
    C2sampledir="$C2sampledir"                  # Comes from the inputfile being read right now.

    # EXPORT FILES
    dirOUT="$HOME/export/$C1exportdir"          # All OUTPUT goes in folders named in col1.
    # ------------                              # ------------
    bam="$HOME/export/$C1exportdir/bam"         # This is where the BAM and SORTED.BAM files go.
    ctab="$HOME/export/$C1exportdir/ctab"       # This is where the CTAB files go.
    fpkm="$HOME/export/$C1exportdir/fpkm"       # This is where the FPKM files go.
    QCraw="$HOME/export/$C1exportdir/QCraw"     # This is where the RAW FastQC output goes.
    QCtrim="$HOME/export/$C1exportdir/QCtrim"   # This is where the TRIMMED FastQC output goes.
    raw="$HOME/export/$C1exportdir/raw"         # This is where the RAW bcl2fastq2 files go.
    sexOUT="$HOME/export/$C1exportdir/Eif2s3y.tsv" # The output file for the sex determination.
    trim="$HOME/export/$C1exportdir/trim"       # This is where the TRIMMED reads go.

    # IMPORT FILES
    adapters="$HOME/import/adapters.sh"         # File with Illumina pair-end adapters to TRIM.
    hisatidx="$HOME/import/HISAT/mm10/genome"   # Basename (no extension) of the reference genome for HISAT.
    inputfile="$HOME/import/input${suffix}.sh"  # CSV with two columns - $dirOUT, $runpath.
    refannot="$HOME/import/GTF/GRCm38.85.gtf"   # Reference gene annotation for StringTie.
    runpath="$HOME/import/NextSeq/$C2sampledir" # All SEQ DATA to be processed is named in col2.

    # FASTQ FILE NAMING CONVENTIONS
    rawNAMES="[!Undetermined]*.fastq.gz"        # Naming convention of the RAW fastq files for QC.
                                                # The '[!Undetermined]*' ignores 'Undetermined' files.
    trimNAMES="[!Undetermined]*R1*.fastq.gz"    # Naming convention for the RAW fastq files being TRIMMED.
                                                # This focuses on the R1 names, because the Trimmomatic
                                                # parameters will change the R1 to R2 later.
                                                # '[1-9]*R1*' focuses only on R1 files starting with a number.

    #┌─────────────────┐
    #│  RUN SCRIPT(S)  ├─ - - - - - - - - - - - - - - - - - - - - - - -
    #└─────────────────┘

    sh ./02_bcl2fastq.sh "$raw" "$runpath"
    sh ./03_fastqc.sh "$raw" "$QCraw" "$rawNAMES"
    sh ./04_trimmomatic.sh "$raw" "$trim" "$adapters" "$trimNAMES"
    sh ./03_fastqc.sh "$trim" "$QCtrim" "$trimNAMES"
    sh ./05_map_align.sh "$trim" "$bam" "$fpkm" "$ctab" "$hisatidx" "$refannot"
    sh ./06_sexing.sh "$bam" "$sexOUT"

    # ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ - - - - - - -

    # LOGFILE
    loopX=$(($SECONDS - $loop0))                # Loop timer ends here.
    sh ./AL_logs.sh "hrs" "$loopX"              # Log the time the loop has run.

done < $inputfile


#┌─────────────────────────────────────────────────────────────────────┐
#│ END SCRIPT                                                          │
#└─────────────────────────────────────────────────────────────────────┘

# LOGFILE
sh ./AL_logs.sh "ple"





