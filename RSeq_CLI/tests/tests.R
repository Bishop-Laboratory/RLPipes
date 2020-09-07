## Test sample-sheets

## Fastqs
# SE fastq (local)
sampleSheet <- data.frame(
  "experiment" = c("RSeq_CLI/tests/NT2_DRIP_head.fastq")
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test1.csv")
# PE fastq (local)
sampleSheet <- data.frame(
  "experiment" = c("RSeq_CLI/tests/qDRIP_R1.fastq+RSeq_CLI/tests/qDRIP_R2.fastq")
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test2.csv")
# PE fastq (local + with control)
sampleSheet <- data.frame(
  "experiment" = c("RSeq_CLI/tests/qDRIP_R1.fastq+RSeq_CLI/tests/qDRIP_R2.fastq"),
  "control" = c("RSeq_CLI/tests/qDRIP_ctr_R1.fastq+RSeq_CLI/tests/qDRIP_ctr_R2.fastq")
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test3.csv")
# SE fastq (local + with control)
sampleSheet <- data.frame(
  "experiment" = c("RSeq_CLI/tests/TC32_DRIP_head.fastq"),
  "control" = c("RSeq_CLI/tests/TC32_DRIP_ctr_head.fastq")
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test4.csv")

## Bams
# SE fastq (local + without control)
sampleSheet <- data.frame(
  "experiment" = c("RSeq_CLI/tests/NT2_DRIP_head.hg38.experiment.bam")
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test5.csv")
# SE fastq (local + with control)
sampleSheet <- data.frame(
  "experiment" = c("RSeq_CLI/tests/TC32_DRIP_head.fastq"),
  "control" = c("RSeq_CLI/tests/TC32_DRIP_ctr_head.fastq")
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test4.csv")






