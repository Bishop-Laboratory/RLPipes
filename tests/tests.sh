#!/bin/bash

testBasicDryDRIP() {
  RSeq -e SRX1025890_TC32_NT_DRIP.hg38.bam -g hg38 -m DRIP -S dryrun=True -o testBasicDryDRIP
}

testBasicDRIP() {
  RSeq -e SRX1025890_TC32_NT_DRIP.hg38.bam -g hg38 -m DRIP -o testBasicDRIP
}

testInputDRIP() {
  RSeq -e SRX1025890_TC32_NT_DRIP.hg38.bam -g hg38 -m DRIP -c SRX1025893_TC32_Input.hg38.bam -o testInputDRIP
}

# Load shUnit2.
. ./shunit2
