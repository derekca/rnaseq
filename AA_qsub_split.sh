#!/bin/sh
# AUTHOR : Derek Caetano-Anolles
# EDITED : 2016.10
# USAGE  : ./AA_qsub_split.sh
# DESCR. : Takes a two-column input file and splits it into individual
#          (one-line) input files that will each be individually qsub'd.
#          Alternatively, running the submission .sh without this will
#          run the pipeline for each row of the input, one after another.


inputfile="$HOME/import/input.sh"            # CSV with two columns - $dirOUT, $runpath.
                                             # For example:
                                             # cfw01 160322_NS500351_0115_AH5KKCAFXX
                                             # cfw02 160525_NS500351_0135_AH2F2YBGXY
                                             # cfw03 160502_NS500351_0128_AHWT7FBGXX
                                             # cfw04 160503_NS500351_0129_AHY5YVBGXX

while read dirOUT runpath                    # Read the 2 columns from the inputfile.
do

    split="$HOME/import/input${dirOUT}.sh"   # Split job inputs go into the same directory.
    echo "$dirOUT" "$runpath" > $split       # This writes the line into the new split file
                                             # and names it input{suffix}.sh

    qsub -v suffix=$dirOUT -V ./00_submit.sh # Submits the script, along with the "suffix"
                                             # environmental variable.

done < $inputfile
wait                                         # Wait for background jobs to complete.






