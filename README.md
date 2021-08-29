# RSeq CLI

**RSeqCLI** is a quality-centric pipeline for R-loop-mapping data analysis.

The primary outputs of the pipeline are:
1. Coverage (.bw) tracks 
2. Peaks (.broadpeak) files
3. Alignment (.bam) files
4. An HTML report that includes:
    - Quality metrics
    - Genomic feature analysis
    - Gene enrichment analysis
    - Comparative analysis with [RMapDB](https://github.com/Bishop-Laboratory/RMapDB)

## Quickstart

The pipeline uses [snakemake](https://snakemake.readthedocs.io/en/stable/) for
executing workflow jobs, and it uses both [RSeqR](https://github.com/Bishop-Laboratory/RSeqR) 
and [RMapDB](https://github.com/Bishop-Laboratory/RMapDB) for downstream analysis. 

### Install

In the future the RSeq CLI will be installed using conda:

```shell
conda install -c bioconda rseq-cli
```

**Until the first official release on conda, you need to build the package via the following:**

1. Clone the repo

```shell
git clone https://github.com/Bishop-Laboratory/RSeqCLI.git
```

2. Checkout the cleanupCLI branch

```shell
cd RSeqCLI/
git checkout cleanupCLI
```

3. Build the environment and install.

```shell
conda install -c conda-forge mamba -y
mamba env create -f mamba-environment.yml --force
conda activate rseq
pip install -e .
```

4. Set your GitHub dev token as a variable in R. Generate a token [here](https://github.com/settings/tokens).

```R
GITHUB_PAT="GH_TOKEN_HERE"
```

5. Install the `RSeqR` package using R:

```R
remotes::install_github('Bishop-Laboratory/RSeqR', auth_token='$GITHUB_PAT', dependencies=FALSE, lib='$RLIBPATH')
```

### Basic Usage

To run RSeqCLI, you will need a `samples.csv` file that details the samples in your dataset. 
Here is an example file provided for testing purposes:

|experiment|control   |
|----------|----------|
|SRX113814 |          |
|SRX1025890|SRX1025893|
|SRX1025899|          |

The basic usage of RSeq follows a three-step process: **build** the workflow,
**check** the workflow, and **run** the workflow.

#### **Build**

Building the workflow generates the **config.json** file which is used to 
configure the RSeq workflow. 

```shell
RSeqCLI build -m DRIP rseq_out/ tests/test_data/samples.csv
```

Output:

```shell
Success! RSeq has been initialized at the specified directory: rseq_out/

Run 'RSeqCLI check rseq_out/' to verify the configuration.
```

#### **Check**

This tests that the run will succeed with the given `config.json` file. 

```shell
RSeqCLI check rseq_out/
```

Output:

```shell
Success! The DAG has been generated successfully. You can view it here: rseq_out/dag.png

Run 'RSeqCLI run rseq_out/' to execute the workflow.
```

It also produces a visualization of the workflow DAG.


#### **Run**

Executes the workflow rules.

```shell
RSeqCLI run rseq_out/
```

If multiple cores are available, they can be specified using the `--threads/-t` option.

```shell
RSeqCLI run -t 30 rseq_out/
```


## Usage Reference

Top-level usage:

```shell

Usage: RSeqCLI [OPTIONS] COMMAND [ARGS]...

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

Usage: RSeqCLI build [OPTIONS] RUN_DIR SAMPLES

  Configure an RSeq workflow.

  RUN_DIR: Directory for RSeq Execution. Will be created if it does not exist.

  SAMPLES: A CSV file with at least one column "experiment" that provides the
  path to either local fastq files, bam files, or public sample accessions
  (SRX or GSM).  Input controls should be in the "control" column.

  If providing paired-end fastq files, enter: "exp_1.fastq~exp_2.fastq".

  Columns may also include "genome" and "mode" columns. These will  override
  the -g, -m, and -n  options.

  "genome" (-g/--genome) is not required if providing public data accessions.



  Example #1: "RSeqCLI build -m DRIP samples.csv"; where "samples.csv" is:

  experiment

  SRX113812

  SRX113813



  Example #2: "RSeqCLI build samples.csv"; where "samples.csv" is:

  experiment,control,genome,mode

  qDRIP_siGL3_1.fq~qDRIP_siGL3_2.fq,,hg38,qDRIP

  DRIPc_3T3.fq,Input_3T3.fq,mm10,DRIPc

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
Usage: RSeqCLI check [OPTIONS] RUN_DIR

  Validate an RSeq workflow.

  RUN_DIR: Directory configured with `RSeqCLI build` and ready for checking
  and execution.

Options:
  -s, --smargs TEXT      Dict of arguments passed to the snakemake python API.
                         Default: "{'use_conda': True}". Read the snakemake
                         API reference for the full list of options.
  -t, --threads INTEGER  Number of threads to use. Default: 1
  --bwamem2              Align with BWA-MEM2 instead of BWA. BWA MEM2 Needs >
                         70GB RAM avaialble to build index, but shows > 3x
                         speed increase. Default: False.
  --macs3                Call peaks using macs3 instead of macs2. Default:
                         True.
  --debug                Run pipeline on subsampled number of reads (for
                         testing).
  --help                 Show this message and exit.
```

### Run


```shell
Usage: RSeqCLI run [OPTIONS] RUN_DIR

  Execute an RSeq workflow.

  RUN_DIR: Directory configured with `RSeqCLI build` and ready for checking
  and execution.

Options:
  -s, --smargs TEXT      Dict of arguments passed to the snakemake python API.
                         Default: "{'use_conda': True}". Read the snakemake
                         API reference for the full list of options.
  -t, --threads INTEGER  Number of threads to use. Default: 1
  --bwamem2              Align with BWA-MEM2 instead of BWA. BWA MEM2 Needs >
                         70GB RAM avaialble to build index, but shows > 3x
                         speed increase. Default: False.
  --macs3                Call peaks using macs3 instead of macs2. Default:
                         True.
  --debug                Run pipeline on subsampled number of reads (for
                         testing).
  --help                 Show this message and exit.

```


## Development notes

### Testing

To run the RSeq pipeline you will need R-loop mapping data in either `fastq`,
`bam`, `bigWig`, or `bedGraph` format. You may also use public data accessions
with RSeq, from SRA, GEO, or BioProject (e.g., `SRX1025899`, `SRX8908682`, `GSM4714836`, etc).

You can get some example data like so:

```shell
cd tests/
bash setup_tests.sh
```

This will download a DRIP-Seq sample and the corresponding input control (downsampled due to size).

To test the full RSeq pipeline:

```shell
bash tests.sh
```

You can also test it directly using the CLI interface:

```shell
# Run a dryrun test on the example data
RSeq -e SRX1025890_TC32_NT_DRIP.hg38.bam -g hg38 -m DRIP -S dryrun=True -o testBasicDryDRIP
```
