# RSeq

An flexible pipeline for R-loop mapping experiments

## Quick start

As a command-line tool:
```
RSeq -m DRIPSeq -x GSM2452072,GSM2452073 -c GSM2668157,GSM2668158
```

As an R-package:
```
library(RSeq)
sampleSheet <- data.frame(
  "Experiment" = c("GSM2452072", "GSM2452073"),
  "Control" = c("GSM2668157", "GSM2668158")
)
RSeq(mode = "DRIPSeq", outDir = ".", groups = sampleSheet)
```

These equivalent commands will pull fastq files from the study `GSE93368` and:
1. Trim, filter, calculate quality metrics on reads
2. Align reads to the genome and remove duplicates
3. Create coverage tracks
4. Call peaks
5. Assess R-loop mapping quality
6. Generate standard visualizations



