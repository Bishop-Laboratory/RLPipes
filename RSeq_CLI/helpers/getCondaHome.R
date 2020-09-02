getCondaHome <- function(CONDA_PREFIX) {
  return(gsub(x = CONDA_PREFIX,
     pattern = "(.+/miniconda3)/.+", replacement = "\\1"))
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]

# Return result
result <- getCondaHome(arg)
cat(result)