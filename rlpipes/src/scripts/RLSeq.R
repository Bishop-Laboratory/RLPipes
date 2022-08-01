# Open log and redirect stderr
con <- file(snakemake@log[[1]], open="wt")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

# Get configs for run
config <- jsonlite::read_json(
  snakemake@params$config
)
peaks <- snakemake@input$peaks
coverage <- snakemake@input$coverage
i <- which(config$experiment == snakemake@wildcards$sample)[1]

# Print configs
print(
  list(
    Coverage = coverage,
    Peaks = peaks,
    genome =  config$genome[[i]]
  )
)

# Creat RLRanges object
rlr <- RLSeq::RLRanges(
  peaks = peaks,
  coverage = coverage,
  genome = config$genome[[i]],
  mode = config$mode[[i]],
  sampleName = config$name[[i]]
)

# Run RLSeq
rlr <- RLSeq::RLSeq(rlr)

# Save output
save(rlr, file = snakemake@output$data)

# Create RLSeq report
tmp <- file.path(snakemake@wildcards$outdir, "tmp", "report", snakemake@wildcards$sample)
print(tmp)
RLSeq::report(rlr, reportPath = snakemake@output$html, intermediates_dir = tmp)
