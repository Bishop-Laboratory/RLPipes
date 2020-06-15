
# Default analysis params
default_params <- list(
  "download_sra" = list(
    "description" = "Downloads SRA file from NCBI if using public data.",
    "tool" = list(
      "name" = "sra-tools",
      "version" = "2.10.7",
      "repo" = "sra-tools"
    ),
    "command" = "prefetch [options] {current_sample} --output-file {fastq_dir}/{current_sample}.sra",
    "options" = list()
  ),
  "sra_to_fastq" = list(
    "description" = "Converst SRA file to FastQ format",
    "tool" = list(
      "name" = "sra-tools",
      "version" = "2.10.7",
      "repo" = "bioconda"
    ),
    "command" = "fastq-dump [options] {fastq_dir}/{current_sample}.sra | ",
    "options" = list()
  ),
  "trim_filter" = list(
    "tool" = list(
      "name" = "fastp",
      "version" = "0.20.1",
      "repo" = "bioconda"
    ),
    "command" = "fastp [options] | ",
    "options" = list(
      "--stdout",
      "--stdin",
      "-w {cores}",
      "-h {fastq_QC_dir}/{current_sample}.fastp.html",
      "if_list" = list(
        "SE" = "",
        "PE" = "--interleaved_in"
      )
    )
  ),
  "align_reads" = list(
    "description" = "Align reads to genome",
    "tool" = list(
      "name" = "bwa",
      "version" = "0.7.17",
      "repo" = "bioconda"
    ),
    "command" = "bwa mem [options] {genome_dir} [input] | ",
    "options" = list(
      "-t {cores}"
    ),
    "input" = list(
      "if_list" = list(
        "from_pipe" = "-",
        "from_file" = "{fastq_input}"
      )
    )
  ),
  "sam_to_bam" = list(
    "description" = "Convert sam to bam file",
    "tool" = list(
      "name" = "samtools",
      "version" = "1.10",
      "repo" = "bioconda"
    ),
    "command" = "samtools view [options] - | ",
    "options" = list(
      "-b",
      "-@ {cores}"
    )
  ),
  "sort_bam" = list(
    "description" = "Sort bam file by position",
    "tool" = list(
      "name" = "samtools",
      "version" = "1.10",
      "repo" = "bioconda"
    ),
    "command" = "samtools sort [options] -",
    "options" = list(
      "-o {bam_dir}/{current_sample}.bam",
      "-@ {cores}"
    )
  ),
  "index_bam" = list(
    "description" = "Index bam file",
    "tool" = list(
      "name" = "samtools",
      "version" = "1.10",
      "repo" = "bioconda"
    ),
    "command" = "samtools index [options] {bam_dir}/{current_sample}.bam",
    "options" = list(
      "-@ {cores}"
    )
  ),
  "split_strands" = list(
    "description" = "If library is strand-specific, split into strand-specific bams",
    "tool" = list(
      "name" = "samtools",
      "version" = "1.10",
      "repo" = "bioconda"
    ),
    "command" = "[options]",
    "options" = list(
      "if_list" = list(
        "non_strand_specific" = "",
        "strand_specific & PE" = list(
          "samtools view -@ {cores} -b -f 128 -F 16 {bam_dir}/{current_sample}.bam > {bam_dir}/{current_sample}.fwd1.bam",
          "samtools view -@ {cores} -b -f 80 {bam_dir}/{current_sample}.bam > {bam_dir}/{current_sample}.fwd2.bam",
          "samtools merge -@ {cores} -f {bam_dir}/{current_sample}.m.bam {bam_dir}/{current_sample}.fwd1.bam {bam_dir}/{current_sample}.fwd2.bam",
          "samtools index -@ {cores} {bam_dir}/{current_sample}.m.bam",
          "samtools view -@ {cores} -b -f 144 {bam_dir}/{current_sample}.bam > {bam_dir}/{current_sample}.rev1.bam",
          "samtools view -@ {cores} -b -f 64 -F 16 {bam_dir}/{current_sample}.bam > {bam_dir}/{current_sample}.rev2.bam",
          "samtools merge -@ {cores} -f {bam_dir}/{current_sample}.p.bam {bam_dir}/{current_sample}.rev1.bam {bam_dir}/{current_sample}.rev2.bam",
          "samtools index -@ {cores} {bam_dir}/{current_sample}.p.bam"
        ),
        "strand_specific & SE" = list(
          "samtools view -@ {cores} -b -F 16 {bam_dir}/{current_sample}.bam > {bam_dir}/{current_sample}.m.bam",
          "samtools view -@ {cores} -b -f 16 {bam_dir}/{current_sample}.bam > {bam_dir}/{current_sample}.p.bam"
        )
      )
    )
  ),
  "alignment_QC" = list(
    "description" = "Check alignment quality.",
    "tool" = list(
      "name" = "qualimap",
      "version" = "2.2.2a",
      "repo" = "bioconda"
    ),
    "command" = "qualimap bamqc -bam {bam_dir}/{current_sample}.bam [options]",
    "options" = list(
      "-nt {cores}",
      "-outdir {bam_QC_dir}/{current_sample}"
    )
  ),
  "call_peaks" = list(
    "description" = "Call peaks from bam files. If strand-specific, make total, plus, and minus peaks.",
    "tool" = list(
      "name" = "macs2",
      "version" = "2.2.7.1",
      "repo" = "bioconda"
    ),
    "command" = "macs2 callpeak [options]",
    "options" = list(
      "-g {effective_genome_size}",
      "--outdir {peaks_dir}/{experiment_name}",
      "if_list" = list(
        "PE" = "-f BAMPE",
        "DRIP" = "--broad",
        "DRIPc" = "--broad"
      ),
      "if_list" = list(
        # Based on https://github.com/PEHGP/drippipline
        "strand_specific" = list("-n {experiment_name}",
                                 "-t {bam_dir}/{experiment_name}.bam",
                                 "-c {bam_dir}/{control_name}.bam",
                                 "&& {CMD_PRIOR}",
                                 "-n {experiment_name}_minus",
                                 "-t {bam_dir}/{experiment_name}.m.bam",
                                 "-c {bam_dir}/{control_name}.m.bam",
                                 "&& {CMD_PRIOR}",
                                 "-n {experiment_name}_plus",
                                 "-t {bam_dir}/{experiment_name}.p.bam",
                                 "-c {bam_dir}/{control_name}.p.bam"),
        "non_strand_specific" = list("-n {experiment_name}",
                                     "-t {bam_dir}/{experiment_name}.bam",
                                     "-c {bam_dir}/{control_name}.bam")
      )
    )
  ),
  "calculate_coverage" = list(
    "description" = "Calculate signal coverage from bam files. If strand-specific, make split bigWig files.",
    "tool" = list(
      "name" = "deeptools",
      "version" = "3.1.3",
      "repo" = "bioconda"
    ),
    "command" = "bamCoverage -b {bam_dir}/{current_sample}.bam [options]",
    "options" = list(
      "-p {cores}",
      "--ignoreForNormalization chrX chrY chrM",
      "--ignoreDuplicates",
      "--minMappingQuality 30",
      "--binSize 20",
      "--effectiveGenomeSize {effective_genome_size}",
      "if_list" = list(
        "strand_specific" = list("--filterRNAstrand forward",
                                 "-o {bw_dir}/{current_sample}.m.bw",
                                 "&& {CMD_PRIOR}",
                                 "--filterRNAstrand reverse",
                                 "-o {bw_dir}/{current_sample}.p.bw"),
        "non_strand_specific" = list("-o {bw_dir}/{current_sample}.bw")
      )
    )
  ),
  # Requires experiment/control pair has finished -- calculate R-loop mapping accuracy (ROC)
  "compute_mapping_accuracy" = list(
    "description" = "Use ROC analysis to compute accuracy of R-loop mapping.",
    "tool" = list(
      "name" = "RSeq",
      "version" = "{RSeq_version}",
      "repo" = "bioconductor"
    ),
    "command" = "{RSeq::compute_mapping_accuracy}"
  ),
  # Require all files in condition finished at this point
  "compute_replicate_congruency" = list(
    "description" = "Use compute replicate congruency using internal function.",
    "tool" = list(
      "name" = "RSeq",
      "version" = "{RSeq_version}",
      "repo" = "bioconductor"
    ),
    "command" = "{RSeq::compute_replicate_congruency}"
  ),
  # Requires all files have finished
  "summarize_bam_group" = list(
    "description" = "Summarize read counts in bins for downstream QC of current sample group.",
    "tool" = list(
      "name" = "deeptools",
      "version" = "3.1.3",
      "repo" = "bioconda"
    ),
    "multiBamSummary bins [options]",
    "options" = list(
      "-b {all_experiment_bams}",
      "-o {result_dir}/bamSummary.npz",
      "-l {all_experiment_labels}",
      "-p {cores}",
      "--ignoreDuplicates",
      "--minMappingQuality 30"
    )
  ),
  "bam_PCA" = list(
    "description" = "calculate PCA for group of bam files."

  ),
  "compute_matrix" = list(
    "tool" = list(
      "name" = "deeptools",
      "version" = "3.1.3",
      "repo" = "bioconda"
    ),
    "command" = ""
  )

)
# usethis::use_data(default_params)






