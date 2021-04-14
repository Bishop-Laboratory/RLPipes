# RSeq CLI

RSeq CLI is a command-line interface to the RSeq best-practices R-looping mapping pipeline. 

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

## Basic Usage

```
RSeq [-e experiment] [-c control] [-m mode] [-g genome] [-n name]

  -e|--experiment        fq/bam/SRA    Experiment(s). Use -e1, -e2 for paired-end fq files.
  -c|--control           fq/bam/SRA    Control(s). Use -c1, -c2 for paired-end fq files.
  -m|--mode              str           Analysis mode ('DRIP', 'RChIP', etc.)
  -g|--genome            str           UCSC genome for samples (e.g., 'hg38')
  -n|--name              str           Sample names(s) for use in output files.
  -s|--sampleSheet       csv           CSV file describing samples (overrides -e/-c/-m/-g/-n).
  -o|--outputDir         dir           Output directory. [default = 'RSeq_out/']
  -G|--genomeDir         dir           Genome directory. [default = '~/.RSeq_genomes']
  -t|--threads           int           Specify number of threads. [default = 1]
  -r|--rseqVars          json          rseqVars.json file (above options are ignored).
  -S|--snakeArgs         str           Keyword arguments to snakemake.api. [default: use_conda=True]
  -v|--version                         Display version info
  -h|--help                            Display detailed usage

Examples (see detailed usage for more):

RSeq -e treated.fastq -c untreated.fastq -m DRIP -g mm10 -n my_experiment -o RSeq_out/ -t 20

RSeq -e1 treated_A.R1.fastq treated_B.R1.fastq -e2 treated_A.R2.fastq treated_B.R2.fastq -m R-ChIP -g dm6

```

### Details

Output from `RSeq -h`:

```
Program: RSeq (R-loop mapping best-practices pipeline)
Version: 0.9.9 (April 14, 2021)
Contact: Henry Miller <millerh1@uthscsa.edu>


Usage: RSeq [-e <experiment> | -s <sampleSheet.csv>] [options]


Input options:

-e|--experiment       STR    Main experimental sample(s) for the workflow. This can be local fastq/bam
                             files or public sequence accessions (SRA, GEO, BioSample, or BioProject).
                             Multiple entries should be separated by a space. If supplying paired-end
                             fastq files, use -e1, -e2 nomenclature instead [see below].

                             If supplying a .bam file, the aligned genome must match the genome specified
                             to RSeq with the -g flag. Public data accessions will be used to automatically
                             download the fastq files. If supplying control [-c, -c1, -c2 flags], the order
                             and number of arguments must match the files/accessions supplied to -e or
                             -e1 and -e2. If supplying .fastq files, the genome must be specified [-g flag].
                             If a sampleSheet is provided, this is ignored.

-e1|--experiment_R1   FQ     If providing paired-end fq files, the first mate. If provided, -e is ignored.

-e2|--experiment_R2   FQ     If providing paired-end fq files, the second mate.

-c|--control          STR    Main control sample(s) for the workflow. See experiment [above] for the
                             relevant nomenclature. Samples provided by -c must match those provided in -e.
                             Typically, appropriate files/accessions will be the genomic inpute or
                             RNaseH1-treated control sample corresponding to the experiment of interest.
                             This will be used for peak calling primarily. It is possible to pass "NA" here
                             if there is an experiment among those provided in -e which does not have a
                             corresponding control sample. This option is overwritten if a "control" column
                             is provided by a sampleSheet.

-c1|--control_R1      FQ     If providing paired-end fq files, the first mate. If provided, -c is ignored.

-c2|--control_R2      FQ     If providing paired-end fq files, the second mate.

-n|--name             STR    Sample names(s) for use in output files. If not supplied, they will be inferred
                             from the filename(s) given or, in the case of public accesions, the SRA
                             metadata. If supplied, must match the order of the experiment and control
                             input. Names should not contain special characters or spaces.

-s|--sampleSheet      CSV    CSV file describing the samples. If supplied, -e, -e1, and -e2 are ignored.
                             Convenient way to process large numbers of samples. Columns can include:

                             experiment    equivalent to passing the -e flag. If using paired-end fastq
                                           files, supply them using "file_one.R1.fastq+file_two.R2.fastq"
                                           nomenclature.
                             control       Same as passing -c flag. If using paired-end fastq files,
                                           supply them using "file_one.R1.fastq+file_two.R2.fastq"
                                           nomenclature. Overwrites -c, -c1, and -c2.
                             mode          equivalent to passing the -m flag. Overwrites -m flag.
                             genome        equivalent to passing the -g flag. Overwrites -g flag.
                             sample_name   equivalent to passing the -n flag. Overwrites -n flag.

-r|--rseqVars         JSON   JSON file describing samples and run conditions. This is generated
                             automatically by the RSeq if -r is not specified [<out_dir>/rseqVars.json].
                             If supplied, all input and pipeline options are ignored.


Pipeline options:

-o|--outputDir        DIR    Output directory. All workflow files will be generated here.

-g|--genome           STR    UCSC genome if supplying fastq or bam files. It will be inferred if using
                             public samples via SRA metadata taxonomy. If genome is not found in the
                             genomeDir [-G flag], it will be downloaded/indexed automatically. All UCSC
                             genomes may be used. If you wish to use a custom genome, submit an issue on
                             GitHub. Custom support will be added in a future release if requested.

-G|--genomeDir        DIR    Directory where RSeq genomes are stored, defaults to [~/.RSeq_genomes].
                             RSeq searches this directory for the requested genome and build if not found.

-t|--threads          INT    Specify number of threads to be used in parallel computations. [default = 1]

-S|--snakeArgs        STR    Keyword arguments to pass to snakemake.api. For the full list of
                             arguments: https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html
                             By default, use_conda=True. This ensures dependency handeling for the pipeline.
                             E.g., RSeq -e testFile.fq -m DRIP -S touch=True report=True containerize=False


Misc. options:

-v|--version                 Display the version info.

-h|--help                    Display detailed usage.

```

