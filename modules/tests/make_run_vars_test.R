library(RSeq)

# Initialize minimal sample sheet with public data
sample_sheet <- data.frame(
  "experiment" = c("SRX2481503", "GSM2326832")
)
run_vars <- makeRSeqDataSet(mode = "DRIP", samples = sample_sheet)




