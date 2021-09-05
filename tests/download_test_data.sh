#!/bin/bash

# Get data
aws s3 sync s3://rlseq-test-data/bam-files/ tests/test_data/ --no-sign-request --no-progress

# Make fastq files
samtools bam2fq tests/test_data/SRX6427717_hg38.bam -1 tests/test_data/SRX6427717_hg38_1.fq -2 tests/test_data/SRX6427717_hg38_2.fq
samtools bam2fq tests/test_data/SRX6427719_hg38.bam -1 tests/test_data/SRX6427719_hg38_R1.fq -2 tests/test_data/SRX6427719_hg38_R2.fq
samtools bam2fq tests/test_data/SRX1025890_hg38.bam -0 tests/test_data/SRX1025890_hg38.fastq
samtools bam2fq tests/test_data/SRX1025893_hg38.bam -0 tests/test_data/SRX1025893_hg38.fastq
samtools bam2fq tests/test_data/SRX1674681_mm10.bam -0 tests/test_data/SRX1674681_mm10.R1.fq
