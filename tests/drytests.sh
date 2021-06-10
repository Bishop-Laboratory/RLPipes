#!/bin/bash

# testBasicDryDRIP() {
#   bash ../bin/RSeq -e SRX1025890_TC32_NT_DRIP.hg38.bam -g hg38 -m DRIP -S dryrun=True -o testBasicDryDRIP
# }
# 
# testBasicDRIP() {
#   bash ../bin/RSeq -e SRX1025890_TC32_NT_DRIP.hg38.bam -g hg38 -m DRIP -o testBasicDRIP -t 40
# }


# testDryInputDRIP() {
#   bash ../bin/RSeq -e SRX1025890_TC32_NT_DRIP.hg38.bam -g hg38 -m DRIP -c SRX1025893_TC32_Input.hg38.bam -o testDryInputDRIP -S dryrun=True
# }
# 
# testInputDRIP() {
#   bash ../bin/RSeq -e SRX1025890_TC32_NT_DRIP.hg38.bam -g hg38 -m DRIP -c SRX1025893_TC32_Input.hg38.bam -o testInputDRIP -t 40
# }



testDryInputDRIPMACS3() {
  bash ../bin/RSeq -e SRX1025890_TC32_NT_DRIP.hg38.bam -g hg38 -m DRIP -c SRX1025893_TC32_Input.hg38.bam -M -o testDryInputDRIPM3 -S dryrun=True
}

testInputDRIPMACS3() {
  bash ../bin/RSeq -e SRX1025890_TC32_NT_DRIP.hg38.bam -g hg38 -m DRIP -c SRX1025893_TC32_Input.hg38.bam -o testInputDRIPM3 -M -t 40
}


# Load shUnit2.
. ./shunit2
