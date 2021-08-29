#!/bin/bash

# Get data
aws s3 sync s3://rseq-testing/bam-files/ tests/test_data/ --no-sign-request --no-progress
