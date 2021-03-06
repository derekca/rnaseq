# RNASEQ PIPELINE

### *Contents*

* [***I. General Overview***](#anchor-1)
  * *A. File List*
  * *B. Pipeline Overview*

* [***II. Running the Pipeline***](#anchor-2)
  * *A. Required Inputs*
  * *B. Generated Export Files*

* [***III. Author Notes***](#anchor-3)
  * *A. Author*
  * *B. Repository*
  * *C. Licenses*
  * *D. Acknowledgements*

- - -





<a id="anchor-1"></a>

## I. General Overview 


### *A. File List*


- `./AA_qsub_split.sh` — Splits an input file into individual jobs for the 00 script. Not necessary in the pipe, but convenient when there are many jobs. Works by taking a two-column input file (input.sh) and splits it into individual (one-line) input files that will each be individually qsub'd. Alternatively, running the submission .sh without this will run the pipeline for each row of the input, one after another.

- `./00_submit.sh` `<suffix>` — Generic script for running any batch job to a queue. Submits the `01` script to the Cluster Nodes for batch processing. This is the starting point for all scripts run on the cluster.

- `./01_pipeline.sh` `<suffix>` — Runs script for project data. All inputs out output directories are specified in this script. the inputs for all subscripts and it will carry over to the others. 

- `./02_bcl2fastq2.sh` `<raw>` `<QCraw>` `<runpath>` `<rawNAMES>` — Converts the (raw) bcl files to fastq, using bcl2fastq.

- `./03_fastqc.sh` `<out>` `<QCout>` `<QCnames>` — Checks quality of RAW or TRIMMED fastq files.

- `./04_trimmomatic.sh` `<raw>` `<trim>` `<QCtrim>` `<adapters>` `<trimNAMES>` — Trims the adapters off of the RAW fastq files, using Trimmomatic.

- `./05_map_align.sh` `<trim>` `<bam>` `<fpkm>` `<ctab>` `<hisatidx>` `<refannot>` — Maps/aligns the trimmed reads with the reference genome. Uses HISAT and StringTie to map/align reads that are found in a "trim" folder.

- `./06_sexing.sh` `<bam>` `<sexOUT>` — Uses a Y-chromosome gene (Eif2s3y) to confirm the sex of the samples. Reads bam files out of a 'bam' folder and exports the number of genes that appear into a spreadsheet.

- `./AL_logs.sh` `<module>` `<logjob>` `<input1>` `<input2>` — This script is in charge of writing the logfile. The logfile is written in a modular way, with different logfile outputs being written based on which module is called. For example, calling "ini" writes the job initiation output at the start of the logfile.


### *B. Pipeline Overview*


This pipeline takes raw NextSeq data, converts the BCL files to FASTQ files, trims the adapters off of those sequences, and then aligns the trimmed sequences onto a reference genome.

An additional step sexes the (mouse) samples by searching for the presence of a mouse Y-chromosome gene, Eif2s3y. However, this particular step can be excluded in non-mouse samples entirely, or modified with a different gene and chromosomal coordinate for the user's specific use. This step serves only as an additional validation of the samples being sequenced, but does not contribute to the alignments or their analysis in any way.

The pipeline requires a number of necessary files from the user (sequence data, reference genomes, etc), along with an input file containing the names of the samples to be run. These files will be described in the following sections.

The pipeline is initiated by submitting the `00` script to the server (using `qsub`). Since this script runs the pipeline on all of the samples presented to it in a sequential manner, this will be slow for a large number of samples. Instead, the `AA` script can be used to split the job up into individual samples before submitting it to the server under the `00` script. The `00` script is only responsible submitting the actual pipeline script, `01`, and can easily be modified to include other functions. The entire pipeline is contained in `01`, along with the paths for all the input and output directories (all of these paths need to be modified by the user before being run). This script then calls seperate scripts for bcl2fastq2 `02` to convert the raw files to fastq format, FastQC `03` to check the quality of the raw files, and Trimmomatic `04` to trim the adapters off of the raw files. The FastQC `03` script is then called again in order to check the quality of the trimmed files, and then finally the trimmed sequences are mapped and aligned `05` to the reference genome using HISAT2 and StringTie. If using *M. musculus* sequences, then a final step sexes the samples `06` by looking for the Y-chromosomal gene *Eif2s3y* as a way to validate the sex of the input.

These scripts are modular, so it is very easy to subsititute one script for another. For example, the script for HISAT2 mapping could be substituted with BWA-MEM using the same inputs.

Throughout this whole process, the logging script `AL` is called and writes into a logfile to quickly keep track of the status of each script, and for debugging purposes. The scripts will log a 'Failure' for a specific script if no file is being generated by that script, however it is important to note that `AL` will still log a 'Success' if the script generates *any* file, even if that file is nonsense. Therefore, it is always reccomended that the user check the outputs to validate that the pipeline is running appropriately for the specific samples.





<a id="anchor-2"></a>

## II. Running the Pipeline


### *A. Required Inputs*


This pipeline requires several files as input from the user before being ready to run, all of which should be downloaded or generated beforehand. All of these variables should be defined by the user in the `01_pipeline.sh` script, easily located in the appropriately annotated section. Descriptions of each file are below.

```
┬
└─▢ $HOME              ▢ home directory
  ├─▢ ../$runpath      ▢ location of Next-seq data
  └─▢ $import          ▢ location of "Import Files"
    ├─▢ adapters.fa    ▢ adapters file
    ├─▢ input.sh       ▢ input file
    ├─▢ hisatidx.fa    ▢ hisat reference index
    └─▢ refannot.gtf   ▢ stringtie gene annotation reference
```

- `$runpath` — This path directs the script to the folder containing all of the folders containing the raw Next-Seq data. Make sure that each folder in this directory contains a spreadsheet called 'SampleSheet.csv' (case-sensitive) with the sample data, as this is read by bcl2fastq2.

- `adapters.fa` — Points to a file containing all the Illumina pair-end adapters to be trimmed from the raw data. Adapters can be found in the 'SampleSheet.csv' which is contained in the raw NextSeq data files, or the equivalent file for your data. This file should be produced on your own, using the following template:
```
>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
```

- `hisatidx` — Points to directory with the reference genome for HISAT2. More info can be found at the [HISAT2 website](https://ccb.jhu.edu/software/hisat2) along with reference indices that can downloaded, if you do not want to generate your own. Regardless of the index's source, the reference genome being pointed at with this variable should be named 'genome.fa'.

- `input.sh` — Points to a `.fa` file with two columns. Column 1 (also called `C1exportdir` in the script comments) contains the user-selected output directory names, which is something useful for the user to identify which samples are being run. Column 2 (also called `C2sampledir` in the script comments) contains the input directory names, which is the name of the folder being described by `$runpath`. These folders can sometimes have human-unfriendly names, which is why including custom names into Column 1 is permitted. The columns need to be able to be read by a `while read` loop, so they must be separated by a single space or a tab. Leave the final line of the list blank.
```
sample01 160322_NS500351_0115_AH5KKCAFXX
sample02 160525_NS500351_0135_AH2F2YBGXY
sample03 160502_NS500351_0128_AHWT7FBGXX
sample04 160503_NS500351_0129_AHY5YVBGXX

```

- `refannot` — Points to the directory for the reference gene annotation file used by StringTie to guide the assembly process. The [StringTie website][refann1] offers suggestions on where to generate the reference annotation file (in GTF or GFF3 format). The reference annotation file can come from any appropriate source, ie. [UCSC,][refann2] [Ensembl,][refann3] [NCBI,][refann4] etc. Ensembl has them under [*Downloads > Download Data via FTP*][refann5].

[refann1]: https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
[refann2]: http://genome.ucsc.edu/
[refann3]: http://www.ensembl.org/index.html
[refann4]: https://www.ncbi.nlm.nih.gov/
[refann5]: http://www.ensembl.org/info/data/ftp/index.html

Additional inputs from the user include naming conventions for files that need to be read by the various scripts. These include:

- `rawNAMES` — Naming convention of the raw fastq files for use by FastQC. This depends on how the Illumina raw data is named (ie. 'Sample1.fastq.gz', 'Sample 2.fastq.gz', etc) and also which files the user is interested in looking at.

- `trimNAMES` — Naming convention for the trimmed fastq files for use by FastQC. 

In addition, the following tools must be downloaded and installed on your server in order to run the scripts in this repository. The version numbers used to generate these scripts are included, but other versions should work unless substantial changes were made to the tool's accepted syntax.

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


### *B. Generated Export Files*


This pipeline generates a directory containing many sets of files as output. The paths where these directories are made is completely up to the user, but as a default they will be contained within a folder (with the sample name), inside of another folder called 'Export', which is inside of the `$HOME` folder. The following structure is used:

```
┬
└─▢ $HOME              ▢ home directory
  ├─▢ $exportdir       ▢ output path is user-specified
  │ ├─▢ bam            ▢ bam  files
  │ ├─▢ ctab           ▢ ctab files
  │ ├─▢ Eif2s3y.tsv    ▢ sex determination spreadsheet
  │ ├─▢ fpkm           ▢ fpkm files
  │ ├─▢ QCraw          ▢ raw  fastq quality check
  │ ├─▢ QCtrim         ▢ trim fastq quality check
  │ ├─▢ raw            ▢ raw  fastq files from bcl
  │ └─▢ trim           ▢ trim fastq files
  ├─▢ $logdir          ▢ location of logfiles
  │ ├─▢ $logerr        ▢ exported sterr output in its folder
  │ └─▢ $logfile       ▢ logfile in its folder
  ├─▢ stderr.log       ▢ most recent standard error file
  └─▢ stdout.log       ▢ most recent standard output file
```

- `bam` : This is where the BAM and SORTED.BAM files go from HISAT.
- `ctab` : This is where the CTAB files go from StringTie.
- `fpkm` : This is where the FPKM files go from StringTie.
- `QCraw` : This is where the raw FastQC output goes, referenced as `QCout` in the FastQC script.
- `QCtrim` : This is where the trimmed FastQC output goes, referenced as `QCout` in the FastQC script.
- `raw` : This is where the raw bcl2fastq2 files go, referenced as `out` in the FastQC script.
- `sexOUT` : The output file for the sex determination.
- `trim` : This is where the trimmed reads go after Trimmomatic, referenced as `out` in the FastQC script.





<a id="anchor-3"></a>

## III. Author Notes


- **Author:** Derek Caetano-Anolles
- **Website:** [derekca.xyz](http://derekca.xyz)
- **Repository:** [github.com/derekca](https://github.com/derekca)
- **Licenses:** Unless otherwise stated, the materials presented in this package are distributed under the [MIT License.](https://opensource.org/licenses/MIT)
- **Acknowledgements:** These materials are based upon work supported by the National Science Foundation Postdoctoral Research Fellowship in Biology under [Grant No. 1523549.](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1523549) Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

- **Edited:** 2016.10.31






