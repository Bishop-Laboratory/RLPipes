#' Correlation test
#' Test the correlation between the current sample other samples around gold-standard genes
correlation_test <- function(input, sample_name, mode, helper_dir) {

	### Bug testing ##
	#input <- "~/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/RSeq_out16/SRX113812_Ntera2_DNA/QC/SRX113812_Ntera2_DNA.hg38.gold_standard_bin_scores.tab"
	#sample_name <- "SRX113812_Ntera2_DNA"
	#mode <- "DRIP"
	#helper_dir <- "~/Bishop.lab/Projects/RSeq/RSeq_CLI/helpers/"
  ############

	public_vals <- paste0(helper_dir, "/data/gold_standard_bw_coverage.rda")
	md_template <- paste0(helper_dir, "/report_template.Rmd")

	output_correlation_picture <- file.path(dirname(input), gsub(basename(input),
																											pattern = "\\.gold_standard_bin_scores\\.tab",
																											replacement = ".correlation.png"))
	output_data <- file.path(dirname(input), gsub(basename(input),
																											pattern = "\\.gold_standard_bin_scores\\.tab",
																											replacement = ".correlation.rda"))
	output_html <- file.path(dirname(input), gsub(basename(input),
																											pattern = "\\.gold_standard_bin_scores\\.tab",
																											replacement = ".correlation.html"))

	# Wrangle the RMapDB bin counts and the new data together
	load(public_vals)
	covInput <- readr::read_tsv(input)
	inputMat <- as.matrix(dplyr::select(covInput, -1, -2, -3))
	rownames(inputMat) <- paste0(covInput[,1, drop = TRUE], "_", covInput[,2, drop = TRUE], "_", covInput[,3, drop = TRUE])
	inputMat <- inputMat[rownames(mat),, drop = FALSE]
	newMat <- cbind(mat, inputMat)
	colnames(newMat)[length(colnames(newMat))] <- sample_name
	colData$Source <- "RMapDB"
	colnames(colData)[1] <- "Mode"
	annoNow <- rbind(colData, data.frame(Mode = mode, row.names = sample_name, Source = "User-supplied"))
	corMat <- cor(newMat)
	
	## TODO: incorporate this into an interactive report
	#relevant_samples <- colnames(corMat)[grep(colnames(corMat), pattern = paste0("^",mode,"_"))]
	#relevant_samples <- relevant_samples[grep(relevant_samples, invert = TRUE, pattern = "RNaseH1")]
	#corMatRel <- corMat[relevant_samples,relevant_samples]
	#avg <- mean(corMatRel[,length(colnames(corMatRel))])
	# Do distance and then downsmaple for static pictures

	# Make static correlation picture
	# From: https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
	paletteLength <- 100
	myColor <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(paletteLength)
	# length(breaks) == length(paletteLength) + 1
	# use floor and ceiling to deal with even/odd length pallettelengths
	myBreaks <- c(seq(min(corMat), 0, length.out=ceiling(paletteLength/2) + 1), 
	              seq(max(corMat)/paletteLength, max(corMat), length.out=floor(paletteLength/2)))
	
	# Plot the heatmap
	pheatmap::pheatmap(corMat, show_rownames = TRUE,
										 show_colnames = TRUE, 
										 color = myColor, breaks = myBreaks,
										 height = 40, width = 40,
										 annotation_colors = list(
										   Source = c("RMapDB" = "grey", "User-supplied" = "firebrick")
										 ),
										 annotation_col = annoNow,
										 filename = output_correlation_picture)
	save(corMat, annoNow, file = output_data)
	# # Make markdown report TODO: Needs to include more info...
	# rmarkdown::render(md_template, params = list(data_file = output_data), output_format = "html_document",
	# 									output_file = output_html)

}

# Parse shell args
arg <- commandArgs(trailingOnly=TRUE)

# Get output JSON
suppressMessages(correlation_test(input = arg[1],
													 sample_name = arg[2],
													 mode = arg[3],
													 helper_dir = arg[4]))
