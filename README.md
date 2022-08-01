# RLPipes
<img src="https://rlbase-data.s3.amazonaws.com/misc/assets/whitebgRLPipes+Logo.png" align="right" alt="logo" width="240" style = "border: none; float: right;">

![Build Status](https://github.com/Bishop-Laboratory/RLPipes/workflows/tests/badge.svg) [![codecov](https://codecov.io/gh/Bishop-Laboratory/RLPipes/branch/main/graph/badge.svg)](https://codecov.io/gh/Bishop-Laboratory/RLPipes) ![Version](https://anaconda.org/bioconda/rlpipes/badges/version.svg) ![license](https://anaconda.org/bioconda/rlpipes/badges/license.svg) ![downloads](https://anaconda.org/bioconda/rlpipes/badges/downloads.svg) 


**RLPipes** is an upstream workflow for R-loop-mapping data. 

The primary outputs of the pipeline are:
1. Coverage (.bw) tracks 
2. Peaks (.broadpeak) files
3. Alignment (.bam) files
4. [RLSeq](https://github.com/Bishop-Laboratory/RLSeq) report (.html and .rda) files

Following RLPipes, the [RLSeq](https://github.com/Bishop-Laboratory/RLSeq) R 
package can be used for more fine-grained downstream analysis.

## Install

The preferred installation method is `mamba` or `conda` (slower):

```shell
mamba create -n rlpipes -c bioconda -c conda-forge rlpipes
conda activate rlpipes
```

### Using `pip`

RLPipes can also be installed with `pip`. However, system dependencies will 
still need to be installed. To accomplish this, do the following:

```shell
git clone https://github.com/Bishop-Laboratory/RLPipes.git
cd RLPipes/
conda install -c conda-forge mamba -y
mamba env create -f rlpipes.yml --force
conda activate rlpipes
python -m pip install -e .
```

## Basic Usage

To run RLPipes, you will need a `samples.csv` file that describes your samples. 
Here is an example file provided for testing purposes:

|experiment|control   |
|----------|----------|
|SRX113814 |          |
|SRX1025890|SRX1025893|
|SRX1025899|          |

The basic usage of RSeq follows a three-step process: **build**, **check** , and **run**.

### **Build**

`RLPipes build` generates a **config.json** file that controls the underlying `snakemake` workflow.

```shell
RLPipes build -m DRIP rlpipes_out/ tests/test_data/samples.csv
```

Output:

```shell
Success! RSeq has been initialized at the specified directory: rlpipes_out/

Run 'RLPipes check rlpipes_out/' to verify the configuration.
```

### **Check**

Verifies that the run will succeed and generates a plot of the workflow jobs. 

```shell
RLPipes check rlpipes_out/
```

Output:

```shell
Success! The DAG has been generated successfully. You can view it here: rlpipes_out/dag.png

Run 'RLPipes run rlpipes_out/' to execute the workflow.
```

### **Run**

Executes the workflow rules.

```shell
RLPipes run rlpipes_out/
```

If multiple cores are available, they can be specified using the `--threads/-t` option.

```shell
RLPipes run -t 30 rlpipes_out/
```

## Usage Reference

Top-level usage:

```shell

Usage: RLPipes [OPTIONS] COMMAND [ARGS]...

  RSeq: An R-loop mapping pipeline with built-in QC.

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  build  Configure an RSeq workflow.
  check  Validate an RSeq workflow.
  run    Execute an RSeq workflow.
  
```

### Build

```shell
Usage: RLPipes build [OPTIONS] RUN_DIR SAMPLES

  Configure an RLPipes workflow.

  RUN_DIR: Directory for RLPipes Execution. Will be created if it does not
  exist.

  SAMPLES: A CSV file with at least one column "experiment" that provides the
  path to either local fastq files, bam files, or public sample accessions
  (SRX or GSM). Input controls should be in the "control" column.

  If providing paired-end fastq files, enter: "exp_1.fastq~exp_2.fastq".

  Columns may also include "genome" and "mode" columns. These will override
  the -g, -m, and -n  options.

  "genome" (-g/--genome) is not required if providing public data accessions.



  Example #1: "RLPipes build -m DRIP outdir/ samples.csv"

  samples.csv:

      experiment

      SRX113812

      SRX113813



  Example #2: "RLPipes build outdir/ samples.csv"

  samples.csv:

      experiment, control, genome, mode

      qDRIP_siGL3_1.fq~qDRIP_siGL3_2.fq, , hg38, qDRIP

      DRIPc_3T3.fq, Input_3T3.fq, mm10, DRIPc



Options:
  -m, --mode TEXT    The type of sequencing (e.g., "DRIP"). The available
                     options are currently: DRIP, DRIPc, qDRIP, sDRIP, ssDRIP,
                     R-ChIP, RR-ChIP, RDIP, S1-DRIP, DRIVE, RNH-CnR, and MapR
  -g, --genome TEXT  UCSC genome for samples (e.g., 'hg38'). Not required if
                     providing public data accessions.
  -n, --name TEXT    Sample names for use in output report. By default,
                     inferred from inputs.
  --help             Show this message and exit.
```

### Check

```shell
Usage: RLPipes check [OPTIONS] RUN_DIR

  Validate an RLPipes workflow.

  RUN_DIR: Directory configured with `RLPipes build` and ready for checking
  and execution.

Options:
  -s, --smargs TEXT      Dict of arguments passed to the snakemake python API.
                         Default: "{'use_conda': True}". Read the snakemake
                         API reference for the full list of options.
  -t, --threads INTEGER  Number of threads to use. Default: 1
  --bwamem2              Align with BWA-MEM2 instead of BWA. BWA MEM2 Needs >
                         70GB RAM avaialble to build index, but shows > 3x
                         speed increase. Default: False.
  --macs2                Call peaks using macs2 instead of macs2
  -G, --groupby TEXT     Column(s) which identify biologically-meaningful
                         grouping(s) of samples (i.e., conditions).  Can be
                         any column name from the `samples` CSV file. If using
                         public data accessions,  it may also include "study".
                         NOTE: If --groupby is set and there R-loop-mapping
                         and expression samples within groups, expression-
                         matched analysis will be run. This can be disabled
                         with the --noexp flag.
                         
                         Example #1: "RSeqCLI build outdir/ samples.csv
                         --groupcols tissue"
                         
                             samples.csv:
                         
                               experiment, mode, tissue
                         
                               GSM1720615, DRIP, NT2
                         
                               GSM1720616, DRIP, NT2
                         
                               GSM1720619, DRIP, K562
                         
                         
                         
                           Example #2: "RSeqCLI build outdir/ samples.csv
                          --groupby tissue"
                         
                             samples.csv:
                         
                               experiment, mode, tissue
                         
                               GSM1720615, DRIP, NT2
                         
                               GSM1720616, DRIP, NT2
                         
                               GSM1720613, DRIPc, NT2
                         
                               GSM1720614, DRIPc, NT2
                         
                               GSM1720622, RNA-seq, NT2
                         
                               GSM1720623, RNA-seq, NT2
                         
  --noexp                If set, no expression-matched analysis will be
                         performed.
  --noreport             If set, RLSeq reports will not be generated.
  --debug                Run pipeline on subsampled number of reads (for
                         testing).
  --tsv                  Obtain config from config.tsv file instead of
                         config.json.
  --noaws                If set, prefetch from SRA tools will be used to 
                         download any public SRA data instead of AWS S3.
  --help                 Show this message and exit.
```

### Run

```shell
Usage: RLPipes run [OPTIONS] RUN_DIR

  Execute an RLPipes workflow.

  RUN_DIR: Directory configured with `RLPipes build` and ready for checking
  and execution.

Options:
  -s, --smargs TEXT      Dict of arguments passed to the snakemake python API.
                         Default: "{'use_conda': True}". Read the snakemake
                         API reference for the full list of options.
  -t, --threads INTEGER  Number of threads to use. Default: 1
  --bwamem2              Align with BWA-MEM2 instead of BWA. BWA MEM2 Needs >
                         70GB RAM avaialble to build index, but shows > 3x
                         speed increase. Default: False.
  --macs2                Call peaks using macs2 instead of macs2
  -G, --groupby TEXT     Column(s) which identify biologically-meaningful
                         grouping(s) of samples (i.e., conditions).  Can be
                         any column name from the `samples` CSV file. If using
                         public data accessions,  it may also include "study".
                         NOTE: If --groupby is set and there R-loop-mapping
                         and expression samples within groups, expression-
                         matched analysis will be run. This can be disabled
                         with the --noexp flag.
                         
                         Example #1: "RSeqCLI build outdir/ samples.csv
                         --groupcols tissue"
                         
                             samples.csv:
                         
                               experiment, mode, tissue
                         
                               GSM1720615, DRIP, NT2
                         
                               GSM1720616, DRIP, NT2
                         
                               GSM1720619, DRIP, K562
                         
                         
                         
                           Example #2: "RSeqCLI build outdir/ samples.csv
                          --groupby tissue"
                         
                             samples.csv:
                         
                               experiment, mode, tissue
                         
                               GSM1720615, DRIP, NT2
                         
                               GSM1720616, DRIP, NT2
                         
                               GSM1720613, DRIPc, NT2
                         
                               GSM1720614, DRIPc, NT2
                         
                               GSM1720622, RNA-seq, NT2
                         
                               GSM1720623, RNA-seq, NT2
                         
  --noexp                If set, no expression-matched analysis will be
                         performed.
  --noreport             If set, RLSeq reports will not be generated.
  --debug                Run pipeline on subsampled number of reads (for
                         testing).
  --tsv                  Obtain config from config.tsv file instead of
                         config.json.
  --help                 Show this message and exit.
```
