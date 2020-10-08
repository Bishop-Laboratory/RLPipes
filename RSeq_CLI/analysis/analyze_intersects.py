import pybedtools
from pybedtools import BedTool

tiles_10kb = BedTool("RSeq_CLI/analysis/10kb_tiles_hg38.bed")
current_bed = BedTool("RSeq_CLI/tests/RSeq_out16/SRX255724_NT2_DRIP-seq_2_rep_2/peaks_macs_unstranded/SRX255724_NT2_DRIP-seq_2_rep_2_hg38_peaks.broadPeak")

tiles_10kb.intersect(current_bed).saveas("RSeq_CLI/analysis/10kb_tiles_hg38_intersected.bed")



