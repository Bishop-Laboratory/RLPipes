#!/bin/bash

## Test local fastq files


RSeq -m DRIP -s RSeq_CLI/tests/sampleSheet_test1.csv -g hg38 -o RSeq_CLI/tests/RSeq_out1 -t 96

RSeq -m qDRIP -s RSeq_CLI/tests/sampleSheet_test2.csv -g hg38 -o RSeq_CLI/tests/RSeq_out2 -t 96

RSeq -m qDRIP -s RSeq_CLI/tests/sampleSheet_test3.csv -g hg38 -o RSeq_CLI/tests/RSeq_out3 -t 96

RSeq -m DRIP -s RSeq_CLI/tests/sampleSheet_test4.csv -g hg38 -o RSeq_CLI/tests/RSeq_out4 -t 96


## Starting from BAM

# SE DRIP
RSeq -m DRIP -s RSeq_CLI/tests/sampleSheet_test5.csv -g hg38 -o RSeq_CLI/tests/RSeq_out5 -t 96

# SE DRIP + control
RSeq -m DRIP -s RSeq_CLI/tests/sampleSheet_test6.csv -g hg38 -o RSeq_CLI/tests/RSeq_out6 -t 96

# PE stranded DRIP
RSeq -m DRIP -s RSeq_CLI/tests/sampleSheet_test7.csv -g hg38 -o RSeq_CLI/tests/RSeq_out7 -t 96

# PE stranded DRIP + control
RSeq -m DRIP -s RSeq_CLI/tests/sampleSheet_test8.csv -g hg38 -o RSeq_CLI/tests/RSeq_out8 -t 96

## Test public data


## RChIP with fastq files


## Starting from bigWigs



## Mismatched bigWig

