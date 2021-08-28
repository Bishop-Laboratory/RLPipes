#!/bin/bash

# Get data
aws s3 sync s3://rseq-testing/bam-files/ test_data/
