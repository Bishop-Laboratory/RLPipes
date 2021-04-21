# Datasets used in tests
These datasets are large enough for proper pipeline testing, but still small enough to fit on GitHub. For fastq files, 
they are gzipped at the highest compression level available using `gzip -9 <some_file>`. Here is the list of datasets:

1. `SRR393964__drip_s96_single_end.fastq.gz` - fastq file containing 100,000 DRIP-Seq single-end reads from an S9.6 sample
   (ie. should be mapping R-loops)
2. `SRR5137197__drip_input_paired_end.fastq` - fastq file containing 100,000 DRIP-Seq paired-end reads from an input sample
   (ie. should NOT be mapping R-loops, negative control used in peak calling)
3. 
