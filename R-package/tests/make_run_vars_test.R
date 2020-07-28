library(RSeq)

# Simplest init
sample_sheet <- data.frame(
  "experiment" = c("SRX2481503", "GSM2326832")
)
run_vars <- initialize_run(mode = "DRIP", samples = sample_sheet, output_json = "output.json")


# Initialize minimal sample sheet with public data
sample_sheet <- data.frame(
  "control" = c("", "", "", ""),
  "experiment" = c("SRX2481503", "GSM2326832", "", "")
)
run_vars <- initialize_run(mode = "DRIP", samples = sample_sheet, output_json = "output.json")


sample_sheet <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/instance/uploads/wduzweiptxfmtavettrtuenyzozylsvffhglffvjdbobhhwokd/samples.csv"
run_vars <- initialize_run(mode='DRIP', samples=sample_sheet, output_json = "output.json")

