# RSeq CLI

**RSeqCLI** is a quality-centric pipeline for R-loop-mapping data analysis.

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

1. Set your GitHub dev token as an environmental variable. Generate a token [here](https://github.com/settings/tokens).

```shell
GITHUB_PAT="GH_TOKEN_HERE"
```

2. Clone the repo

```shell
git clone https://github.com/Bishop-Laboratory/RSeq.git
```

3. Build the conda recipe in a new environment (`rseq`) and install

```shell
cd RSeq/
conda install -c conda-forge mamba -y
mamba env create -f mamba-environment.yml --force
conda activate rseq
conda mambabuild -c conda-forge -c bioconda bioconda-recipe-testing/ |& tee build.log
BINARY_PATH=$(grep -i "TEST END" build.log | awk '{ print $3 }')
mamba remove rseq  # Remove previous version
conda install $BINARY_PATH
RLIBPATH=$(Rscript -e "cat(.libPaths()[1])")
Rscript -e "remotes::install_github('Bishop-Laboratory/RSeqR', auth_token='$GITHUB_PAT', dependencies=FALSE, lib='$RLIBPATH')"
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

#### Build

Building the workflow generates the **config.json** file which is used to 
configure the RSeq workflow. 

```shell
RSeqCLI build -o rseq_out/ -m DRIP tests/test_data/samples.csv
```

Output:

```shell
Success! The config file is located here: rseq_out/config.json

Run 'RSeqCLI check rseq_out/config.json' to verify the configuration.
```

#### Check

This tests that the run will succeed with the given `config.json` file. 

```shell
RSeqCLI check rseq_out/config.json
```

Output:

```shell
Success! The DAG has been generated successfully. You can view it here: rseq_out/dag.png

Run 'RSeqCLI run rseq_out/config.json' to execute the workflow.
```

It also produces a visualization of the workflow DAG.


#### Run

Executes the workflow rules.

```shell
RSeqCLI run rseq_out/config.json
```

If multiple cores are available, they can be specified using the `--threads/-t` option.

```shell
RSeqCLI run -t 30 rseq_out/config.json
```


## Usage Reference

User manual goes here.

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
