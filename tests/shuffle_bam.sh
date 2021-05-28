#!/bin/bash

bamToBed -i "/home/millerh1/PycharmProjects/RSeq/rseq_out_testfiles/SRX1025890_TC32_NT_DRIP.hg38.bam" -bed12 -cigar | shuffleBed -i - -g hg38.chrom.sizes | bedToBam -bed12 -g hg38.chrom.sizes | samtools sort > "/home/millerh1/PycharmProjects/RSeq/rseq_out_testfiles/SRX1025890_TC32_NT_DRIP.shuffle.hg38.bam"
