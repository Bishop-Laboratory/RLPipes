# RSeq CLI

RSeq CLI is a command-line interface to the RSeq best-practices R-looping mapping pipeline. It's good software!

## Quickstart

RSeq CLI can be installed using conda:

```bazaar
conda install -c bioconda rseq-cli
```

To run the RSeq pipeline you will need R-loop mapping data in either `fastq`,
`bam`, `bigWig`, or `bedGraph` format. You may also use public data accessions
with RSeq, from SRA, GEO, or BioProject (e.g., `SRX8908682`, `GSM4714836`, etc).

To run the basic RSeq workflow:

```
RSeq -m DRIP -s sampleSheet.csv -o RSeq_out/ -g ~/.RSeq_genomes -t 98
```

This will run the full RSeq pipeline to:
1. Deduplicate and condense the fastq files
2. Map the reads to the genome.
3. Call R-loop peaks and signal tracks
4. Calculate the quality score of the resulting maps.
5. Return an analysis report in HTML format.

\* Upon first use, RSeq genomes will not be available and may take some time to be generated.

## Usage

```
RSeq [-m mode] [-s sampleSheet] [-o outputDir] [-g genomeDir] [-t threads]

General options:

  -m|--mode              mode      Choose analysis mode ('DRIP', 'RChIP', etc.) See full list in detailed usage.
  -s|--sampleSheet       file      CSV file describing samples. See details.
  -i|--inputDir          dir       Directory containing fastq/bam files. Not required for SRA samples.
  -o|--outputDir         dir       Project directory. [default = 'RSeq_out/']
  -g|--genomeDir         dir       Genome directory. [default = '~/.RSeq_genomes']
  -t|--threads           int       Specify number of threads. [default = 1]
  --noFastp                        Do not use fastp to perform adapter trimming, filtering, and fastq QC.
  --noDedupe                       Do not use clumpify.sh to perform read deduplication.
  --keepTmp                        Do not delete intermediate files once processing is finished.
  --returnBamsOnly                 Returns only bams, coverage tracks, and normalized signal tracks.
  --noMerge                        Do not attempt to merge technical replicates from SRA.
  --help                           Display detailed usage

```

### Details

Here are some details.
### How to setup a virtual environment with conda
The conda package and environment manager is included in all versions of
[Anaconda](https://docs.conda.io/projects/conda/en/latest/glossary.html#anaconda-glossary),
[Miniconda](https://docs.conda.io/projects/conda/en/latest/glossary.html#miniconda-glossary),
and [Anaconda Repository](https://docs.continuum.io/anaconda-repository/).
Here we will go with the lightweight miniconda option.
+ Download the latest [miniconda installer](https://docs.conda.io/en/latest/miniconda.html) for your operating system and run the installer script.
For example, as of Apr 2021 the latest for MacOS is `Miniconda3-py39_4.9.2-MacOSX-x86_64.sh`. Open a terminal, navigate to the directory where the installer
 got downloaded, and run the installer with bash:
```
bash Miniconda3-py39_4.9.2-MacOSX-x86_64.sh
```
You will be prompted to read the user licence agreement which you can skip to
the end by typing "q". There you will be asked to accept the licence terms.
Go ahead and initialize Miniconda3 when prompted. You will need to close your
terminal and open a new one for the changes to take effect. Go ahead and do so.
By default the conda base environment will be activated whenever you open a new
terminal. You will know the conda base environment is active because you will
see a `(base)` prefix in front of your `user@computer` address.
I prefer to deactivate this default behavior because I have other ongoing
projects and I don't want to interfere with their python settings. To deactivate
the default behavior type:
```
conda config --set auto_activate_base false
```
Next time you open an new terminal you should not see the `(base)` prefix.
You can activate the conda base environment manually with the command
```
conda activate base
```

+ Double check conda is updated.
```
conda update -n base conda
```
+ Create a virtual environment for RSeq with the specified requirements.
```
conda create --file requirements.txt -c conda-forge -n Rseq
```
When prompted new packages will be installed, procceed by pressing `y`.
+ Activate `Rseq` virtual environment
```
conda activate Rseq
```
+ Update virtual environment (when necessary)
```
conda env update -n Rseq -f requirements.txt
```
+ To deactivate environment, use
```
conda deactivate
```
+ To remove an environment, use
```
conda env remove -n Rseq
```

### Continuous Integration
+ Check for code style and run tests by running the integration script.
```
bash integration.sh
```
