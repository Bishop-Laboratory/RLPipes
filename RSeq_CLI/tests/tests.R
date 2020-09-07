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
# TODO: This will produce an error if the filename contains .hg38.experiment.bam
sampleSheet <- data.frame(
  "experiment" = c("RSeq_CLI/tests/TC32_DRIP_head.bam"),
  "control" = c("RSeq_CLI/tests/TC32_DRIP_ctr_head.bam")
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test6.csv")
# PE fastq (local)
# TODO: This will produce an error if the filename contains .hg38.experiment.bam
sampleSheet <- data.frame(
  "experiment" = c("RSeq_CLI/tests/qDRIP.bam")
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test7.csv")
# SE fastq (local + with control)
# TODO: This will produce an error if the filename contains .hg38.experiment.bam
sampleSheet <- data.frame(
  "experiment" = c("RSeq_CLI/tests/qDRIP.bam"),
  "control" = c("RSeq_CLI/tests/qDRIP_ctr.bam")
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test8.csv")

# SE fastq (local + with control)
# TODO: This will produce an error if the filename contains .hg38.experiment.bam
sampleSheet <- data.frame(
  "experiment" = "RSeq_CLI/tests/qDRIP.bam",
  "control" = "RSeq_CLI/tests/qDRIP_ctr.bam"
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test8.csv")


## Public
# SE fastq (without control)
sampleSheet <- data.frame(
  "experiment" = "SRR2019278"
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test9.csv")
# SE fastq (with control)
sampleSheet <- data.frame(
  "experiment" = "SRR2019278",
  "control" = "SRR2019281"
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test10.csv")
# PE fastq (without control)
sampleSheet <- data.frame(
  "experiment" = "SRR10916579"
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test11.csv")
# SE fastq (local + with control)
sampleSheet <- data.frame(
  "experiment" = "SRR10916579",
  "control" = "SRR10916580"
)
write.csv(sampleSheet, "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/sampleSheet_test12.csv")

