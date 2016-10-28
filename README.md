# I. GENERAL OVERVIEW

## *File List*

- *00_submit.sh* — Generic script for running any batch job to a queue. Submits the `01` script to the Cluster Nodes for batch processing. This is the starting point for all scripts run on the cluster.
- *01_pipeline.sh* — Runs script for project data. All inputs out output directories are specified in this script. the inputs for all subscripts and it will carry over to the others. 
- *02_bcl2fastq2.sh* — Converts the (raw) bcl files to fastq, using bcl2fastq.
- *03_fastqc.sh* — Checks quality of RAW or TRIMMED fastq files.
- *04_trimmomatic.sh* — Trims the adapters off of the RAW fastq files, using Trimmomatic.
- *05_map_align.sh* — Maps/aligns the trimmed reads with the reference genome. Uses HISAT and StringTie to map/align reads that are found in a "trim" folder.
- *06_sexing.sh* — Uses a Y-chromosome gene (Eif2s3y) to confirm the sex of the samples. Reads bam files out of a 'bam' folder and exports the number of genes that appear into a spreadsheet.
- *AA_qsub_split.sh* — Splits an input file into individual jobs for the 00 script. Not necessary in the pipe, but convenient when there are many jobs. Works by taking a two-column input file (input.sh) and splits it into individual (one-line) input files that will each be individually qsub'd. Alternatively, running the submission .sh without this will run the pipeline for each row of the input, one after another.
- *AL_logs.sh* — This script is in charge of writing the logfile. The logfile is written in a modular way, with different logfile outputs being written based on which module is called. For example, calling "ini" writes the job initiation output at the start of the logfile.

## *Software Used*

The following tools must be downloaded and installed on your server in order to run the scripts in this repository. The version numbers used to generate these scripts are included, but other versions should work unless substantial changes were made to the tool's accepted syntax.

| Tool                 | Version        | Summary               |
|:-------------------- |:-------------- |:--------------------- |
| [bcl2fastq2][BC]     | 2.17.1.14      | Convert bcl to fastq. |
| [FastQC][QC]         | 0.11.4         | Check fastq quality.  |
| [HISAT2][H2]         | 2.0.1-beta     | Map reads.            |
| [SAMtools][SA]       | 0.1.19-44428cd | Sort alignments.      |
| [StringTie][ST]      | 1.2.1          | Sequence assembly.    |
| [Trimmomatic][TR]    | 0.35           | Trim reads.           |

[BC]: http://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html
[QC]: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[H2]: https://ccb.jhu.edu/software/hisat2/index.shtml
[SA]: http://samtools.sourceforge.net/
[ST]: https://ccb.jhu.edu/software/stringtie/
[TR]: http://www.usadellab.org/cms/?page=trimmomatic



# II. SCRIPT PIPELINE SUMMARY

## *A. Script Usage*

- *./00_submit.sh* `<suffix>`

- *./01_pipeline.sh* `<suffix>`

- *./02_bcl2fastq2.sh* `<raw>` `<QCraw>` `<runpath>` `<rawNAMES>`

- *./03_fastqc.sh* `<out>` `<QCout>` `<QCnames>`

- *./04_trimmomatic.sh* `<raw>` `<trim>` `<QCtrim>` `<adapters>` `<trimNAMES>`

- *./05_map_align.sh* `<trim>` `<bam>` `<fpkm>` `<ctab>` `<hisatidx>` `<refannot>`

- *./06_sexing.sh* `<bam>` `<sexOUT>`

- *./AA_qsub_split.sh*

- *./AL_logs.sh* `<module>` `<logjob>` `<input1>` `<input2>`

```
READ FROM INPUT.SH, DON'T EDIT
C1exportdir="$C1exportdir"                  # Comes from the inputfile being read right now.
C2sampledir="$C2sampledir"                  # Comes from the inputfile being read right now.

EXPORT FILES

IMPORT FILES

FASTQ FILE NAMING CONVENTIONS
rawNAMES     Naming convention of the RAW fastq files for QC.
              The '[!Undetermined]*' ignores 'Undetermined' files.
trimNAMES  Naming convention for the RAW fastq files being TRIMMED.
            This focuses on the R1 names, because the Trimmomatic
              parameters will change the R1 to R2 later.
                '[1-9]*R1*' focuses only on R1 files starting with a number.
```



# III. FILES AND VARIABLES

## *A. Import Files*

The following files are necessary to run these scripts, and should be assembled beforehand. Descriptions of each file are below.

```
DIRECTORY              DESCRIPTION
---------------------------------------------------------------
┬
└─▣ $HOME              ▣ home directory
  ├─▣ ../$runpath      ▣ location of seq-data
  └─▣ $import          ▣ location of "Import Files"
    ├─▣ adapters.fa    ▣ adapters file
    ├─▣ input.sh       ▣ input file
    ├─▣ hisatidx.fa    ▣ hisat reference index
    └─▣ refannot.gtf   ▣ stringtie gene annotation reference

```


- `adapters.fa` — Points to a file containing all the adapters to be trimmed. Adapters can be found in the `SampleSheet.csv` which is contained in the raw NextSeq data files, or the equivalent file for your data. This file should be produced on your own, using the following template:
```
>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
```

- `hisatidx` — Points to directory with the reference genome for HISAT2. More info can be found at the [HISAT2 website](https://ccb.jhu.edu/software/hisat2) along with reference indices that can downloaded, if you do not want to generate your own. Regardless of the index's source, the reference genome being pointed at with this variable should be named `genome.fa`.

- `input.sh` — Points to .fa file with the output dir names, and the input dir names. Col. 1 are the output folder names, Col. 2 are the input folder names. The columns need to be able to be read by a "while read" loop, so they must be separated by a single space or a tab. Leave the final line blank or it ignores the last entry.
```
cfw01 160322_NS500351_0115_AH5KKCAFXX
cfw02 160525_NS500351_0135_AH2F2YBGXY
cfw03 160502_NS500351_0128_AHWT7FBGXX
cfw04 160503_NS500351_0129_AHY5YVBGXX

```

- `refannot` — Points to directory for the reference gene annotation file used by StringTie.

- `runpath` — All SEQ DATA to be processed is named in $C2sampledir.


adapters          # File with Illumina pair-end adapters to TRIM.
hisatidx    # Basename (no extension) of the reference genome for HISAT.
inputfile   # CSV with two columns - $dirOUT, $runpath.
refannot    # Reference gene annotation for StringTie.
runpath     # All SEQ DATA to be processed is named in col2.





















## *B. Directories/files that are created by these scripts*


```
DIRECTORY            USR¹  DESCRIPTION
---------------------------------------
┬
└─▣ $HOME              ▣    home directory
  │
  ├─▢ $exportdir       ▢    output path is user-specified
  │ ├─▢ bam            ▢    bam  files
  │ ├─▢ ctab           ▢    ctab files
  │ ├─▢ Eif2s3y.tsv    ▢    sex determination spreadsheet
  │ ├─▢ fpkm           ▢    fpkm files
  │ ├─▢ QCraw          ▢    raw  fastq quality check
  │ ├─▢ QCtrim         ▢    trim fastq quality check
  │ ├─▢ raw            ▢    raw  fastq files from bcl
  │ └─▢ trim           ▢    trim fastq files
  │
  ├─▢ $logdir          ▢    location of logfiles
  │ └─▢ $logfile       ▢    logfile in its folder
  │
  ├─▢ stderr.log       ▢    standard error file
  └─▢ stdout.log       ▢    standard output file

```

DIRECTORIES FROM INPUT.SH
-----------------------------
C1exportdir       Points to the output folders that get made to store all the data.
C2sampledir       Points to the directory where the sequence data to be analyzed is located.

EXPORT FILES
-----------------------------
bam               This is where the BAM and SORTED.BAM files go.
ctab              This is where the CTAB files go.
fpkm              This is where the FPKM files go.
QCraw             This is where the RAW FastQC output goes.
QCtrim            This is where the TRIMMED FastQC output goes.
raw               This is where the RAW bcl2fastq2 files go.
sexOUT            The output file for the sex determination.
trim              This is where the TRIMMED reads go.

FASTQ FILE NAMING CONVENTIONS
-----------------------------
rawNAMES          Naming convention of the RAW fastq files for QC.
trimNAMES         Naming convention for the RAW fastq files being TRIMMED.











# IV. SOFTWARE NOTES

## *A. HISAT*

```
hisat2 -x <index>     Run HISAT2 using this reference genome index. Can be downloaded
                      from https://ccb.jhu.edu/software/hisat2
-1                    Comma-separated list of files containing mate 1s.
-2                    Comma-separated list of files containing mate 2s
--rna-strandedness    Specify strand-specific information: the default
                      is unstranded. For single-end reads [F] means a
                      read corresponds to a transcript, and [R] is rev compl.
                      For paired-end reads, use [FR] or [RF].
-dta                  Report alignments tailored for transcript assemblers, ie StringTie.
-p <int>              Launch NTHREADS parallel build threads (default: 1).

samtools sort         Sort alignments by leftmost coordinates, or by read name when the
                      [-n] option is used.
samtools view         With no options or regions specified, prints all alignments in the
                      specified input alignment file (in SAM, BAM, or CRAM format) to
                      standard output in SAM format (with no header). Using [-S] was required
                      (in previous versions) if input was in SAM format, but now the correct
                      format is automatically detected. The [-u] option outputs uncompressed
                      BAM. This option saves time spent on compression/decompression and is
                      thus preferred when the output is piped to another samtools command.
```

## *B. StringTie*

```
stringtie <input.bam> Run StringTie using the input *.bam file(s)
-o <out.gtf>          Sets output GTF name where StringTie writes assembled transcripts.
-p <int>              Specify number of processing threads (CPUs). Default is 1.
-G <ref_ann.gff>    	Use the reference annotation file (in GTF or GFF3 format) to guide
                      the assembly process. The output will include expressed reference
                      transcripts as well as any novel transcripts that are assembled.
                      This option is required by options -B, -b, -e, -C (see below).
                      Can come from UCSC, Ensembl, NCBI, etc. Ensembl has them under
                      "Downloads" > "Download Data via FTP".
-b <path>             Enables the output of *.ctab files for Ballgown, but these files will
                      be created in the provided directory <path> instead of the directory
                      specified by the -o option. Note: adding the -e option is recommended
                      with the -B/-b options, unless novel transcripts are still wanted in
                      the StringTie GTF output.
-e                    Limits the processing of read alignments to only estimate and output
                      the assembled transcripts matching the reference transcripts given
                      with the -G option (requires -G, recommended for -B/-b). Boosts speed.
```

## *C. Trimmomatic*


```
java -jar <path to trimmomatic.jar> \                 Run JAR filepath
          [PE] \                                      Paired End Mode (as opposed to SE)
          [-threads <threads>] \                      indicates # threads to on multi-core computers.
          [-phred33 | -phred64] \                     Specifies the base quality encoding.
          [-trimlog<logFile>] \                       Creates a log of all read trimmings with:
                                                      - read name
                                                      - surviving sequence length
                                                      - loc of 1st surviving base (amt trimmed from start)
                                                      - loc of last surviving base in original read
                                                      - amount trimmed from the end
          [-basein <in> | <inR1><inR2>] \             Input .gz files from the RAW reads.
          [-baseout <out> | <R1p><R1up><R2p><R2up> \  The 4 outputs.
          [...]                                       Any of the additional options, below.

          TRIMMOMATIC OPTIONS:
          ILLUMINACLIP:   Cut adapter and other illumina-specific sequences from the read.
          SLIDINGWINDOW:  Performs a sliding window trimming approach. It starts scanning
                          at the 5‟ end and clips the read once the average quality within
                          within the window falls below a threshold.
          MAXINFO:        An adaptive quality trimmer which balances read length and error
                          rate to maximise the value of each read.
          LEADING:        Cut bases off the start of a read, if below a threshold quality
          TRAILING:       Cut bases off the end of a read, if below a threshold quality
          CROP:           Cut the read to a specified length by removing bases from the end
          HEADCROP:       Cut the specified number of bases from the start of the read
          MINLEN:         Drop the read if it is below a specified length
          AVGQUAL:        Drop the read if the average quality is below the specified level
          TOPHRED33:      Convert quality scores to Phred-33
          TOPHRED64:      Convert quality scores to Phred-64
```







# V. AUTHOR NOTES

- **Author:** Derek Caetano-Anolles
- **Repository:** [github.com/derekca](https://github.com/derekca)
- **Licenses:** Unless otherwise stated, the materials presented in this package are distributed under the [MIT License.](https://opensource.org/licenses/MIT)
- **Acknowledgements:** These materials are based upon work supported by the National Science Foundation Postdoctoral Research Fellowship in Biology under [Grant No. 1523549.](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1523549) Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
- **Edited:** 2016.10.27






