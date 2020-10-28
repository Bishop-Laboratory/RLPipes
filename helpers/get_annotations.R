# Helper script to parse the Homer genome annotation file

args <- commandArgs(trailingOnly = TRUE)
genome <- args[1]
anno_file <- args[2]

anno_folder <- normalizePath(dirname(anno_file))
outdir <- paste0(anno_folder, "/hg38_annos/")

library(rtracklayer)
csfile <- paste0(anno_folder, "/", "chrom.sizes")
download.file(paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/bigZips/", genome, ".chrom.sizes"),
              destfile = csfile)
cs <- read.table(file = csfile, 
                 stringsAsFactors = FALSE)
annos <- readr::read_tsv(anno_file, col_names = FALSE)

cats <- unique(annos$X6)
cats <- cats[! cats %in% c("LTR?", "DNA?", "Unknown", "RC?", "SINE?", "")]
good_cats <- data.frame(
  "acronyms" = cats,
  "meanings" = c("Intergenic", "Simple repeat", "Satellite", "Promoter", "Pseudo gene", "Intron", "TTS",
                 "LINE", "LTR", "SINE", "DNA Tranpsoson", "CpG", "ncRNA-misc", "Low complexity",
                 "Exon", "snRNA", "Retrotransposon", "3'UTR", "5'UTR", "tRNA", "Rolling Circle",
                 "srpRNA", "rRNA", "scRNA", "miRNA", "RNA-misc", "snoRNA"),
  "type" = c("gene", "rep", "rep", "gen", "RNA", "gene", "gene", 
             "rep", "rep", "rep", "rep", "gene", "RNA", "rep",
             "gene", "RNA", "rep", "gene", "gene", "RNA", "rep",
             "RNA", "RNA", "RNA", "RNA", "RNA", "RNA"), stringsAsFactors = FALSE
)

annoshd <- annos
annoshd <- annoshd[,c(2, 3, 4, 1, 5, 6)]
colnames(annoshd) <- c("seqnames", "start", "end", "name", "strand", "class")

for (class_now in cats) {
  print(class_now)
  grNow <- annoshd[annoshd$class %in% class_now,]
  grNow <- GenomicRanges::makeGRangesFromDataFrame(grNow, keep.extra.columns = TRUE)
  seqlevels(grNow, pruning.mode="coarse") <- cs$V1
  genome(grNow) <- genome
  seqlengths(grNow) <- cs$V2
  grNow <- trim(grNow)
  grNow <- sortSeqlevels(grNow)
  grNow <- sort(grNow)
  dir.create(outdir, showWarnings = FALSE)
  export(grNow, con = paste0(outdir, class_now, ".unsorted.bed"))
  system(paste0("bedtools sort -g ", csfile,
                " -i ", paste0(outdir, class_now, ".unsorted.bed") , " > ", paste0(outdir, class_now, ".bed")))
  file.remove(paste0(outdir, class_now, ".unsorted.bed"))
}





