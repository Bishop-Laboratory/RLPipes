library(affy)
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE45544&format=file", destfile = "GSE45544_RAW.tar")
untar("GSE45544_RAW.tar", exdir = "GSE45544_RAW")
AFF <- affy::ReadAffy(filenames = list.files("GSE45544_RAW", full.names = TRUE))
exp_norm <- rma(AFF)
data <- exp_norm@assayData$exprs

library(GEOquery)
res <- GEOquery::getGEO("GSE45544")
colData <- pData(exp)
colData$group <- ""
colData$group[colData$geo_accession %in% paste0("GSM110833", c(4, 6, 7, 8, 9))] <- "Relapse"
colData$group[colData$geo_accession %in% paste0("GSM153466", c(4:8))] <- "Primary"
rmInd <- which(colData$group != "")
colData <- colData[rmInd,]


good_cols <- grep(colnames(data), pattern = paste0(colData$geo_accession, collapse = "|"))
data <- data[,good_cols]
colnames(data) <- colData$geo_accession

library(ggpubr)
library(tidyr)
to_plot <- as.data.frame(data) %>% gather()
ggboxplot(to_plot, x = "key", y = "value", outlier.shape = NA) + rotate_x_text() + scale_y_continuous(limits = c(0, 10))


library(limma)
design <- model.matrix(~0 + group, data = colData)
fit <- lmFit(data, design = design)
contr <- makeContrasts(groupRelapse - groupPrimary, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, n = Inf)
which(abs(top.table$logFC) > 1 & top.table$P.Value < .05)

annot <- read.table("GPL6244-17930.txt", sep = "\t", fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
dfl <- lapply(annot$ID, FUN = function(id_now) {
  ens <- gsub(paste0(" ", annot$gene_assignment[annot$ID == id_now]), pattern = ".+(ENST[0-9]+).+", replacement = "\\1")
  data.frame(id_now, ens)
})
dfl <- data.table::rbindlist(dfl)
dfl <- dfl[grep(dfl$ens, pattern = "^ENST[0-9]+$"),]

library(EnsDb.Hsapiens.v86)
columns(EnsDb.Hsapiens.v86)
tx2sym <- select(EnsDb.Hsapiens.v86, keys = as.character(dfl$ens), keytype = "TXID",  columns = c("SYMBOL", "TXID"))
dfl2 <- merge(x = dfl, y = tx2sym, by.x = "ens", by.y = "TXID")

top.table$id_now <- as.numeric(rownames(top.table))
top.table <- merge(x = dfl2, y = top.table, by = "id_now")

top.table$pval <- -log10(top.table$P.Value)
library(EnhancedVolcano)
EnhancedVolcano(toptable = top.table, lab = "", x = "logFC", y = "P.Value", pCutoff = .05, FCcutoff = 1)
EnhancedVolcano(toptable = top.table, lab = "", x = "logFC", y = "adj.P.Val", pCutoff = .05, FCcutoff = 1)





