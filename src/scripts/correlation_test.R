#' Correlation test
#' Test the correlation between the current sample other samples around gold-standard genes
correlation_test <- function(input, sample_name, mode, helper_dir, output_image, output_data) {
	# Wrangle the RMapDB bin counts and the new data together
	load(paste0(helper_dir, "/data/gold_standard_bw_coverage.rda"))
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
							   filename = output_image)
	save(corMat, annoNow, file = output_data)
}

# Parse shell args
arg <- commandArgs(trailingOnly=TRUE)

# Get output JSON
suppressMessages(correlation_test(input = arg[1],
								  sample_name = arg[2],
								  mode = arg[3],
								  helper_dir = arg[4],
                                  output_image = arg[5],
                                  output_data = arg[6]))
