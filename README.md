# CONTENTS

01  Software used

02  Script pipeline summary

  a  Script pipeline
  b  Script descriptions

03  Files and variables

  a  Directories that are made or used
  b  Passed variable descriptions
  c  Example import file formats

04  Software notes

  a  HISAT
  b  StringTie
  c  Trimmomatic

# I. GENERAL OVERVIEW

## *File List*

- ***00_submit.sh*** — sd
- ***01_pipeline.sh*** — sd
- ***02_bcl2fastq2.sh*** — sd
- ***03_fastqc.sh*** — sd
- ***04_trimmomatic.sh*** — sd
- ***05_map_align.sh*** — sd
- ***06_sexing.sh*** — sd
- ***AA_qsub_split.sh*** — sd
- ***AL_logs.sh*** — sd

## *Software Used*

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

`./00_submit.sh` `<suffix>`

`./01_pipeline.sh` `<suffix>`

`./02_bcl2fastq2.sh` `<raw>` `<QCraw>` `<runpath>` `<rawNAMES>`

`./03_fastqc.sh` `<raw>` `<QCraw>` `<rawNAMES>`

`./04_trimmomatic.sh` `<raw>` `<trim>` `<QCtrim>` `<adapters>` `<trimNAMES>`

`./03_fastqc.sh` `<trim>` `<QCtrim>` `<trimNAMES>`

`./05_map_align.sh` `<trim>` `<bam>` `<fpkm>` `<ctab>` `<hisatidx>` `<refannot>`

`./06_sexing.sh` `<bam>` `<sexOUT>`

`./AA_qsub_split.sh`

`./AL_logs.sh` `<module>` `<logjob>` `<input1>` `<input2>`



```
     DIRECTORY                 PASSED VARIABLES
┬
├─▢ AA_qsub_split.sh
│ │
│ └─▢ 00_submit.sh            suffix
│   │
│   └─▢ 01_pipeline.sh        suffix
│     │
│     ├─▢ 02_bcl2fastq2.sh    raw    QCraw   runpath   rawNAMES
│     │ │
│     │ └─▢ 03_fastqc.sh      raw    QCraw   rawNAMES
│     │
│     ├─▢ 04_trimmomatic.sh   raw    trim    QCtrim    adapters  trimNAMES
│     │ │
│     │ └─▢ 03_fastqc.sh      trim   QCtrim  trimNAMES
│     │
│     ├─▢ 05_map_align.sh     trim   bam     fpkm   ctab   hisatidx   refannot
│     │
│     └─▢ 06_sexing.sh        bam    sexOUT
│
└─▢ AL_logs.sh                module logjob input1     input2
```



## *B. Script Descriptions*

```
─┤AA├─  Splits an input file into individual jobs for the 00 script. Not
 │  │      necessary in the pipe, but convenient when there are many jobs.
─┤00├─  Submits 01 to the Cluster Nodes for batch processing. This is the
 │  │      starting point for all scripts run on the cluster.
─┤01├─  Inputs and outputs for the pipeline, runs the pipe-scripts.
 │  │
─┤02├─  Converts the (raw) bcl files to fastq.
 │  │
─┤03├─  Checks quality of raw (or trimmed) fastq files.
 │  │
─┤04├─  Trims the adapters off of the RAW fastq.
 │  │
─┤05├─  Maps/aligns the trimmed reads with the reference genome.
 │  │
─┤06├─  Uses a Y-chromosome gene to confirm the sex of the samples.
 │  │
─┤AL├─  Updates the logfile when a subscript calls for it.
 │  │
```


# III. FILES AND VARIABLES

## *A. Directories That Are Made Or Used*


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
  ├─▣ $import          ▣    location of "Import Files"
  │ ├─▣ adapters.fa    ▣    adapters file
  │ ├─▣ input.sh       ▣    input file
  │ ├─▣ hisatidx.fa    ▣    hisat reference index
  │ └─▣ refannot.gtf   ▣    stringtie gene annotation reference
  │
  ├─▢ $logdir          ▢    location of logfiles
  │ └─▢ $logfile       ▢    logfile in its folder
  │
  ├─▣ ../$runpath      ▣    location of seq-data
  │
  ├─▢ stderr.log       ▢    standard error file
  └─▢ stdout.log       ▢    standard output file

---------------------------------------
¹Items w/ a filled box are user-required. Others are script-created.
```

## *A. Passed Variable Descriptions


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

IMPORT FILES
-----------------------------
adapters          Points to .fa file with all the adapters to be trimmed.
hisatidx          Points to directory with the reference genome for HISAT2.
                  More info at https://ccb.jhu.edu/software/hisat2
inputfile         Points to .fa file with the output dir names, and the input dir names.
refannot          Points to directory for the reference gene annotation file used by StringTie.
runpath           All SEQ DATA to be processed is named in $C2sampledir.


## C. Example Import File Formats


adapters.sh   Note¹: Adapters can be found in the SampleSheet.csv files of the raw data.

```
>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
```

inputfile.sh  Note¹: Always leave the final line blank or it ignores the last entry.
              Note²: Col. 1 are the output folder names, Col. 2 are the input folder names.
              Note³: The columns need to be able to be read by a "while read" loop,
                     so they must be separated by a single space or a tab.

```
cfw01 160322_NS500351_0115_AH5KKCAFXX
cfw02 160525_NS500351_0135_AH2F2YBGXY
cfw03 160502_NS500351_0128_AHWT7FBGXX
cfw04 160503_NS500351_0129_AHY5YVBGXX
```


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






