########################################################################################################################
##############################################   Parse inputs    #######################################################
########################################################################################################################

import math

# Global Configs
helpers_dir=config['helpers_dir'][0]
genome_home_dir=config['genome_home_dir'][0]
cores=config['threads'][0]
outdir=config['outdir'][0]
outdir=outdir.strip("/")
bwa_mem2=config['bwa_mem2'][0]

# Sample info
mode=config['mode']
paired_end=config['paired_end']
strand_specific=config['strand_specific']
homer_anno_available=config['homer_anno_available']
genome=config['genome']
effective_genome_size=config['effective_genome_size']
full_genome_length=config['full_genome_length']
effective_genome_fraction = [x/y for x, y in zip(effective_genome_size, full_genome_length)]
sample_type=config['file_type']
sample=config['sample_name']
experiments=[exp.split(",") for exp in config['experiment']]
controls=[ctr if ctr != "None" else None for ctr in config['control']]

# Generate the output file names
report_output = expand("{outdir}/RSeq_report/{sample}_{genome}__RSeq_Report.html",zip,
            sample=sample,outdir=[outdir for i in range(len(sample))],genome=genome)

# For testing the workflow on GitHub
test_output="my_file.txt"
test=False
if test:
    report_output = test_output

    rule make_test:
        output: test_output
        conda: helpers_dir + "/envs/bwa.yaml"
        shell: "echo Hello world! > {output}"

# Select bwa type
# BWA MEM2 is still in development and has a particularly problematic habit of over-zealous RAM usage
# https://github.com/bwa-mem2/bwa-mem2/issues/118
# User needs to have the option to choose classic BWA
if bwa_mem2:
    bwa_cmd="bwa-mem2"
    bwa_location="bwa_mem2_index"
else:
    bwa_cmd="bwa"
    bwa_location="bwa_index"

### Helper functions for pipeline ###
def find_replicates(wildcards):
    srr_list = [experiments[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    replicates = expand("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.fastq",
           outdir=wildcards.outdir, srr_acc=srr_list, sample=wildcards.sample)
    return replicates

def pe_test_fastp(wildcards):
    pe = [paired_end[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    if pe:
        res="--interleaved_in "
    else:
        res=""
    return res

def pe_test_samblaster(wildcards):
    pe = [paired_end[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    if not pe:
        res="--ignoreUnmated "
    else:
        res=""
    return res

def pe_test_bwa(wildcards):
    pe = [paired_end[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    if pe:
        res="-p "
    else:
        res=""
    return res

def input_test_callpeak(wildcards):
    input = [controls[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    # print(input)
    if input is not None:
        return_dict = {
            'treatment': wildcards.outdir + "/bam/" + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam",
            'control': wildcards.outdir + "/bam/" + input + "/" + input + "." + wildcards.genome + ".bam",
            'index_treatment': wildcards.outdir + "/bam/" + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam.bai",
            'index_control': wildcards.outdir + "/bam/" + input + "/" + input + "." + wildcards.genome + ".bam.bai"
        }
    else:
        return_dict = {
            'treatment': wildcards.outdir + "/bam/" + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam",
            'control': wildcards.outdir + "/bam/" + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam",
            'index_treatment': wildcards.outdir + "/bam/" + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam.bai",
            'index_control': wildcards.outdir + "/bam/" + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam.bai"
        }

    return return_dict

#######################################################################################################################
##############################################   Main pipeline    ######################################################
########################################################################################################################

rule output:
    input:
        report_output

rule prepare_report:
    input:
        rlfs_enrichment="{outdir}/RLFS_analysis/{sample}/{sample}_{genome}__rlfs_enrichment.rda",
        bam_stats="{outdir}/bam_stats/{sample}/{sample}_{genome}__bam_stats.txt",
        homer_annotations="{outdir}/homer_annotations/{sample}_{genome}__feature_overlaps.txt",
        correlation_analysis="{outdir}/correlation_analysis/{sample}/{sample}_{genome}__correlation_analysis.rda",
        read_qc_data="{outdir}/QC/fastq/json/{sample}.{genome}.json"
    output:
        html="{outdir}/RSeq_report/{sample}_{genome}__RSeq_Report.html"
    conda: helpers_dir + "/envs/prepare_report.yaml"
    params:
        helpers_dir=helpers_dir,
        configs = "{outdir}/config.json",
    log: "{outdir}/logs/prepare_report/{sample}_{genome}__prepare_report.log"
    shell:
     """
     (Rscript {params.helpers_dir}/prepare_report.R {params.helpers_dir} {params.configs} {output} {input}) &> {log}
     """

rule correlation_analysis:
    input: "{outdir}/bin_scores/{sample}/{sample}_{genome}__bin_scores.tab"
    output:
        image="{outdir}/correlation_analysis/{sample}/{sample}_{genome}__correlation_analysis.png",
        data="{outdir}/correlation_analysis/{sample}/{sample}_{genome}__correlation_analysis.rda"
    params:
        helpers_dir = helpers_dir,
        mode = mode
    conda: helpers_dir + "/envs/correlation_analysis.yaml"
    log: "{outdir}/logs/correlation_analysis/{sample}_{genome}__correlation_analysis.log"
    shell:
         """
         (Rscript {params.helpers_dir}/correlation_test.R {input} {wildcards.sample} {params.mode} {params.helpers_dir}) &> {log}
         """

rule calculate_bin_scores:
    input: "{outdir}/coverage/{sample}/{sample}_{genome}__coverage.bw"
    params:
        helpers_dir = helpers_dir
    output:
        npz="{outdir}/bin_scores/{sample}/{sample}_{genome}__bin_scores.npz",
        tab="{outdir}/bin_scores/{sample}/{sample}_{genome}__bin_scores.tab"
    threads: cores
    conda: helpers_dir + "/envs/deeptools.yaml"
    log: "{outdir}/logs/bin_scores/{sample}_{genome}__get_correlation_bin_sthreads.log"
    shell:
        """
        (multiBigwigSummary BED-file --BED {params.helpers_dir}/data/correlation_genes_100kb.{wildcards.genome}.1kbwindow.bed \
        -o {output.npz} -b {input} --outRawCounts {output.tab} -p {threads}) &> {log}
        """

rule calculate_coverage:
    input:
        bam="{outdir}/bam/{sample}/{sample}.{genome}.bam",
        bai="{outdir}/bam/{sample}/{sample}.{genome}.bam.bai"
    output: "{outdir}/coverage/{sample}/{sample}_{genome}__coverage.bw"
    threads: cores
    conda: helpers_dir + "/envs/deeptools.yaml"
    log: "{outdir}/logs/{sample}_{genome}__deeptools_coverage.log"
    params:
        extra="--ignoreForNormalization chrX chrY chrM --minMappingQuality" \
              + " 20 --binSize 10 --effectiveGenomeSize " + str(effective_genome_size)
    shell: """
        (bamCoverage -b {input.bam} -p {threads} {params.extra} -o {output}) &> {log}
    """

rule assign_genome_annotations:
    # Calculates the percentage of called peaks overlapping with repeat regions
    input:
        gene_anno=genome_home_dir + "/{genome}/homer_anno.txt",
        peaks="{outdir}/peaks/{sample}/{sample}_{genome}__compiled_peaks.bed"
    output:
        stats_out="{outdir}/homer_annotations/{sample}_{genome}__feature_overlaps.txt"
    log: "{outdir}/logs/homer_annotations/{sample}_{genome}__homer_annotations.log"
    conda: helpers_dir + "/envs/homer.yaml"
    shell: """
        (assignGenomeAnnotation {input.peaks} {input.gene_anno} > {output.stats_out}) &> {log}
    """

rule download_homer_anno:
    output:
        gene_anno=genome_home_dir + "/{genome}/homer_anno.txt"
    params:
        out_dir=genome_home_dir + "/{genome}/homer_anno",
    log: outdir + "/logs/download_homer_anno/{genome}__download_homer.log"
    shell: """
        (wget -O {params.out_dir}.zip http://homer.ucsd.edu/homer/data/genomes/{wildcards.genome}.v6.4.zip
        unzip -d {params.out_dir} {params.out_dir}.zip
        mv {params.out_dir}/data/genomes/{wildcards.genome}/{wildcards.genome}.full.annotation {output}
        rm -rf {params.out_dir} && rm -rf {params.out_dir}.zip) &> {log}
    """

rule bam_stats:
    input: "{outdir}/bam/{sample}/{sample}.{genome}.bam"
    output: "{outdir}/bam_stats/{sample}/{sample}_{genome}__bam_stats.txt"
    threads: 4
    conda: helpers_dir + "/envs/samtools.yaml"
    log: "{outdir}/logs/bam_stats/{sample}_{genome}__bam_stats.log"
    shell: "(samtools flagstat -@ {threads} {input} > {output}) &> {log}"

rule rlfs_enrichment:
    input:
         peaks="{outdir}/peaks/{sample}/{sample}_{genome}__compiled_peaks.bed",
         rlfs=genome_home_dir + "/{genome}/rloop_predictions/{genome}.rlfs.bed"
    output: "{outdir}/RLFS_analysis/{sample}/{sample}_{genome}__rlfs_enrichment.rda"
    params:
        helpers_dir=helpers_dir
    threads: 5
    conda: helpers_dir + "/envs/rlfs_enrichment.yaml"
    log: "{outdir}/logs/rlfs_enrichment/{sample}_{genome}__rlfs_enrichment.log"
    shell: """
     (Rscript {params.helpers_dir}/rlfs_perm_test.R {threads} {wildcards.genome} {input.peaks} {input.rlfs} {output}) &> {log}
     """

rule download_rlfs_annotations:
    output:
        bed=genome_home_dir + "/{genome}/rloop_predictions/{genome}.rlfs.bed"
    shell: """
    wget -O {output} https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/{wildcards.genome}.rlfs.bed
    """

rule compile_peaks:
    input:
        macs2="{outdir}/peaks/{sample}/macs2/{sample}_{genome}__peaks_macs2.broadPeak",
        epic2="{outdir}/peaks/{sample}/epic2/{sample}_{genome}__peaks_macs2.bed"
    output:
        peaks="{outdir}/peaks/{sample}/{sample}_{genome}__compiled_peaks.bed"
    params:
        helpers_dir=helpers_dir
    log: "{outdir}/logs/compile_peaks/{sample}_{genome}__compile_peaks.log"
    shell: """
    (
    Rscript {params.helpers_dir}/compile_peaks.R {wildcards.sample} {input.macs2} {input.epic2} {output.peaks}
    ) &> {log}
    """

rule macs:
    input: unpack(input_test_callpeak)
    output: "{outdir}/peaks/{sample}/macs2/{sample}_{genome}__peaks_macs2.broadPeak"
    log: "{outdir}/logs/macs2/{sample}_{genome}__macs2.log"
    threads: 1
    conda: helpers_dir + "/envs/macs.yaml"
    shell: """
        (
        if [ {input.control} == {input.treatment} ]; then
            echo "No Control file detected -- running MACS2 without a control"
            macs2 callpeak -t {input.treatment} -n called_peaks/{wildcards.sample}
        else
            echo "Control file detected -- running MACS2 with control"
            macs2 callpeak -t {input.treatment} -c {input.control} -n called_peaks/{wildcards.sample}
        fi
        ) &> {log}
    """

rule epic:
    input:
        #info="{outdir}/info/{sample}.{genome}.experiment_predictd.txt",
        unpack(input_test_callpeak)
    output:
        "{outdir}/peaks/{sample}/epic2/{sample}_{genome}__peaks_macs2.bed"
    log: "{outdir}/logs/epic2/{sample}_{genome}__epic2_callpeaks.log"
    conda: helpers_dir + "/envs/epic.yaml"
    shell: """
        if [ {input.control} == {input.treatment} ]; then
            echo "No Control file detected -- running MACS2 without a control"
            epic2 -t {input.treatment} -gn {wildcards.genome} -o {output}
        else
            echo "Control file detected -- running MACS2 with control"
            macs2 callpeak -t {input.treatment} -c {input.control} -gn {wildcards.genome} -o {output}
        fi
    """

# TODO: Try BWA MEM2 (Has stdin option and 1.5-3X speed increase)
rule bwa_mem:
    input:
        bwa_index_done=[
              genome_home_dir + "/{genome}/bwa_index/{genome}.ann",
              genome_home_dir + "/{genome}/bwa_index/{genome}.pac",
              genome_home_dir + "/{genome}/bwa_index/{genome}.amb"],
        reads="{outdir}/tmp/fastqs_trimmed/{sample}.{genome}__trimmed.fastq"
    output:
        bam="{outdir}/bam/{sample}/{sample}.{genome}.bam",
        bai="{outdir}/bam/{sample}/{sample}.{genome}.bam.bai"
    conda: helpers_dir + "/envs/bwa_mem.yaml"
    priority: 15
    log: "{outdir}/logs/bwa/{sample}_{genome}__bwa_mem.log"
    params:
        index=genome_home_dir + "/{genome}/bwa_index/{genome}",
        bwa_extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        bwa_interleaved=pe_test_bwa,
        samblaster_extra=pe_test_samblaster,
        samtools_sort_extra="-O BAM",
        bwa_cmd=bwa_cmd
    threads: 10
    shell: """
        ({params.bwa_cmd} mem -t {threads} {params.bwa_extra} {params.bwa_interleaved}{params.index} {input.reads} | \
        samblaster {params.samblaster_extra}| \
        samtools view -q 10 -b -@ {threads} - | \
        samtools sort {params.samtools_sort_extra} -@ {threads} -o {output.bam} - && \
        samtools index -@ {threads} {output.bam}) &> {log}
     """

rule bwa_index:
    input:
        genome_home_dir + "/{genome}/{genome}.fa"
    output:
        genome_home_dir + "/{genome}/bwa_index/{genome}.ann",
        genome_home_dir + "/{genome}/bwa_index/{genome}.pac",
        genome_home_dir + "/{genome}/bwa_index/{genome}.amb"
    params:
        prefix=genome_home_dir + "/{genome}/" + bwa_location + "/{genome}",
        bwa_cmd=bwa_cmd
    conda: helpers_dir + "/envs/bwa.yaml"
    log:
        genome_home_dir + "/{genome}/" + bwa_location + "/{genome}_bwa_index.log"
    shell:"""
        ({params.bwa_cmd} index -p {params.prefix} {input}) &> {log}
    """

rule download_fasta:
    output:
        genome_home_dir + "/{genome}/{genome}.fa"
    params:
          prefix=genome_home_dir + "/{genome}/bwa_index/{genome}",
    log: outdir + "/logs/download_fasta/{genome}__download_fasta.log"
    shell: """
        (mkdir -p {params.prefix}
        wget -O {output}.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz
        gunzip {output}.gz && echo "Indexing BAM at {params.prefix} ... This may take some time ...") &> {log}
    """

# TODO: Pipe this part
rule fastp:
    input:
        sample="{outdir}/tmp/fastqs_merged/{sample}.{genome}__merged.fastq"
    output:
        trimmed=temp("{outdir}/tmp/fastqs_trimmed/{sample}.{genome}__trimmed.fastq"),
        html="{outdir}/QC/fastq/html/{sample}.{genome}.html",
        json="{outdir}/QC/fastq/json/{sample}.{genome}.json"
    conda: helpers_dir + "/envs/fastp.yaml"
    log: "{outdir}/logs/fastp/{sample}.{genome}__fastp_pe.log"
    priority: 10
    params:
        extra=pe_test_fastp
    threads: 4
    shell: """
    (fastp -i {input} --stdout {params.extra}-w {threads} -h {output.html} -j {output.json} > {output} ) &> {log}
    """

rule merge_replicate_reads:
    input: find_replicates
    output: temp("{outdir}/tmp/fastqs_merged/{sample}.{genome}__merged.fastq")
    params:
        #debug=" | head -n 40000 -",  # TODO: Why does this throw an error when uncommented?
        debug=""
    log: "{outdir}/logs/merge_fastq/{sample}_{genome}__merge_fastq.log"
    shell: """
    (cat {input}{params.debug} > {output}) &> {log}
    """

# This should always produced interleaved fq even if paired end. I don't like this solution, but it should hold.
# reformat.sh (from BBTools) will interleave the paired-end files
rule sra_to_fastq:
    input: "{outdir}/tmp/sras/{sample}/{srr_acc}/{srr_acc}.sra"
    output: temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.fastq")
    conda: helpers_dir + "/envs/sratools.yaml"
    threads: 1
    log: "{outdir}/logs/sra_to_fastq/{sample}_{srr_acc}__sra_to_fastq_pe.log"
    params:
        output_directory="{outdir}/tmp/sras/{sample}/",
        fqdump="--skip-technical --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --split-3 ",
        debug="",
        #debug=" -X 60000"
    shell: """(
    cd {params.output_directory}
    fastq-dump{params.debug} {params.fqdump}-O ../../fastqs_raw/{wildcards.sample}/ {wildcards.srr_acc}
    cd ../../fastqs_raw/{wildcards.sample}/
    if test -f {wildcards.srr_acc}_2.fastq; then
        echo "Paired end -- interleaving"
        reformat.sh in1={wildcards.srr_acc}_1.fastq in2={wildcards.srr_acc}_2.fastq out={wildcards.srr_acc}.fastq overwrite=true
        rm {wildcards.srr_acc}_1.fastq && rm {wildcards.srr_acc}_2.fastq
    else
        echo "Single end -- finished!"
    fi
    ) &> {log}
    """

# TODO: Figure out how to use the pipes
# TODO: Probably need to specify this version of prefetch and/or find alternative to it...
# TODO: Retry if fails due to network error
rule download_sra:
    output: temp("{outdir}/tmp/sras/{sample}/{srr_acc}/{srr_acc}.sra")
    conda: helpers_dir + "/envs/sratools.yaml"
    log: "{outdir}/logs/download_sra/{sample}__{srr_acc}__download_sra.log"
    params:
        output_directory = "{outdir}/tmp/sras/{sample}/"
    threads: math.ceil(cores * .2)
    shell: """
            (
            cd {params.output_directory}
            prefetch {wildcards.srr_acc} -f yes
            ) &> {log}
            """
