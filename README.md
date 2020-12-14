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



