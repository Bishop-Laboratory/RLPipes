# Getting started with conservation analysis

The raw data for this analysis is here: https://uthscsa.box.com/s/l5ehtq8rd2e7u0ctdl4x3ev3qqgl5rcx

That contains all the peak files. The current R code in `conservation_analysis.R` works if that folder I linked is unzipped in the `conservation_analysis/data` directory. 

For `conservation_analysis.R`, the current strategy is:

1. Make windows across the genome (10kb wide and staggered by 2.5kb).
2. Overlap all the peaks with these windows and count the number that overlap. 
3. Determine, based on this, which windows are 'highly-conserved' or not. 

Questions still to answer:
1. Which are the 'highly-conserved' R-loops?
2. How do they vary with respect to mode of analysis? (e.g., DRIP vs R-ChIP?)
3. What genomic features or genes do they correspond to? 
4. Are they conserved between species?

