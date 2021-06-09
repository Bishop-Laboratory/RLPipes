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
report_html = expand("{outdir}/RSeq_report/{sample}_{genome}__RSeq_Report.html",zip,
            sample=sample, outdir=[outdir for i in range(len(sample))],genome=genome)
report_data = expand("{outdir}/RSeq_report/{sample}_{genome}__RSeq_Report.rda",zip,
            sample=sample, outdir=[outdir for i in range(len(sample))],genome=genome)

# For testing the workflow
debug=True

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
def find_sra_replicates(wildcards):
    srr_list = [experiments[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    replicates = expand("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.fastq",
           outdir=wildcards.outdir, srr_acc=srr_list, sample=wildcards.sample)
    return replicates

def find_fq(wildcards):
    fq = [experiments[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    return fq

def check_type_fq(wildcards):
    file_type = [sample_type[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    if file_type == "fastq":
        fq = [experiments[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0][0]
        if fq[-3:] == ".gz":
            # CASE: Fastq GZ available
            return "{outdir}/tmp/fastqs_gunzip/{sample}.{genome}__gunzip.fastq"
        else:
            # CASE: normal Fastq
            return "{outdir}/tmp/fastqs_cp/{sample}.{genome}__cp.fastq"
    else:
        # CASE: Public data accession
        return "{outdir}/tmp/fastqs_merged/{sample}.{genome}__merged.fastq"
    

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

def get_effective_genome_size(wildcards):
    return [effective_genome_size[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]

def get_mode(wildcards):
    return [mode[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]

def get_control(wildcards):
    return [controls[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]

def get_sample_bam(wildcards):
    return [experiments[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]

def input_test_callpeak(wildcards):
    input = [controls[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    st_now = [sample_type[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]

    # This is a hack that changes the expected bam location to switch between the bwamem rule and bam wrangle rule
    if st_now != "bam":
        bam_type = "/bam/"
    else:
        bam_type = "/wrangled_bam/"

    # Set the expected samples for peak calling
    if input is not None:
        if st_now != "public":
            # If not public, input will be a file. Need input sample name instead.
            input = [sample[idx] for idx, element in enumerate(experiments) if element[0] == input][0]
        return_dict = {
            'treatment': wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam",
            'control': wildcards.outdir + bam_type + input + "/" + input + "." + wildcards.genome + ".bam",
            'index_treatment': wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam.bai",
            'index_control': wildcards.outdir + bam_type + input + "/" + input + "." + wildcards.genome + ".bam.bai"
        }
    else:
        return_dict = {
            'treatment': wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam",
            'control': wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam",
            'index_treatment': wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam.bai",
            'index_control': wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam.bai"
        }

    return return_dict

def get_report_inputs(wildcards):
    return_dict = {
        'rlfs_enrichment': wildcards.outdir + '/RLFS_analysis/' + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.genome + "__rlfs_enrichment.rda",
        'homer_annotations': wildcards.outdir + '/homer_annotations/' + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.genome + "__feature_overlaps.txt",
        'correlation_analysis': wildcards.outdir + '/correlation_analysis/' + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.genome + "__correlation_analysis.rda",
        'peak_compilation_data': wildcards.outdir + '/peaks/' + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.genome + "__compiled_peaks.rda",
        'peaks': wildcards.outdir + '/peaks/' + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.genome + "__compiled_peaks.bed"
    }
    st_now = [sample_type[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    if st_now in ['fastq', 'public', 'bam']:
        return_dict['bam_stats'] = wildcards.outdir + '/bam_stats/' + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.genome + "__bam_stats.txt"
    elif st_now in ['fastq', 'public']:
        return_dict['read_qc_data'] = wildcards.outdir + '/QC/fastq/json/' + wildcards.sample + "." + wildcards.genome + ".json"

    return return_dict

def choose_bam_type(wildcards):
    st_now = [sample_type[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    if st_now != "bam":
        bam_type = "/bam/"
    else:
        bam_type = "/wrangled_bam/"
    return wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam"


#######################################################################################################################
##############################################   Main pipeline    ######################################################
########################################################################################################################

rule output:
    input:
        html=report_html,
        data=report_data

rule prepare_report:
    input: unpack(get_report_inputs)
    output:
        html="{outdir}/RSeq_report/{sample}_{genome}__RSeq_Report.html",
        data="{outdir}/RSeq_report/{sample}_{genome}__RSeq_Report.rda"
    conda: helpers_dir + "/envs/prepare_report.yaml"
    params:
        helpers_dir=helpers_dir,
        configs = "{outdir}/config.json",
    log: "{outdir}/logs/prepare_report/{sample}_{genome}__prepare_report.log"
    shell:
     """
     (Rscript {params.helpers_dir}/scripts/prepare_report.R {params.helpers_dir} {params.configs} {wildcards.sample} \
     {output.html} {output.data} {input}) &> {log}
     """

rule correlation_analysis:
    input: "{outdir}/bin_scores/{sample}/{sample}_{genome}__bin_scores.tab"
    output:
        image="{outdir}/correlation_analysis/{sample}/{sample}_{genome}__correlation_analysis.png",
        data="{outdir}/correlation_analysis/{sample}/{sample}_{genome}__correlation_analysis.rda"
    params:
        helpers_dir = helpers_dir,
        mode = get_mode
    conda: helpers_dir + "/envs/correlation_analysis.yaml"
    log: "{outdir}/logs/correlation_analysis/{sample}_{genome}__correlation_analysis.log"
    shell:
         """
         (Rscript {params.helpers_dir}/scripts/correlation_test.R {input} {wildcards.sample} {params.mode} {params.helpers_dir} {output.image} {output.data}) &> {log}
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
    input: choose_bam_type
    output: "{outdir}/coverage/{sample}/{sample}_{genome}__coverage.bw"
    threads: cores
    conda: helpers_dir + "/envs/deeptools.yaml"
    log: "{outdir}/logs/coverage/{sample}_{genome}__deeptools_coverage.log"
    params:
        effective_genome_size=get_effective_genome_size,
        extra="--ignoreForNormalization chrX chrY chrM --minMappingQuality" \
              + " 20 --binSize 10 --effectiveGenomeSize "
    shell: """
        (bamCoverage -b {input} -p {threads} {params.extra}{params.effective_genome_size} -o {output}) &> {log}
    """

rule assign_genome_annotations:
    # Calculates the percentage of called peaks overlapping with repeat regions
    input:
        gene_anno=genome_home_dir + "/{genome}/homer_anno.txt",
        peaks="{outdir}/peaks/{sample}/{sample}_{genome}__compiled_peaks.bed"
    output:
        stats_out="{outdir}/homer_annotations/{sample}/{sample}_{genome}__feature_overlaps.txt"
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
    input: choose_bam_type
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
     (Rscript {params.helpers_dir}/scripts/rlfs_perm_test.R {threads} {wildcards.genome} {input.peaks} {input.rlfs} {output}) &> {log}
     """

rule download_rlfs_annotations:
    output:
        bed=genome_home_dir + "/{genome}/rloop_predictions/{genome}.rlfs.bed"
    shell: """
    wget -O {output} https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/{wildcards.genome}.rlfs.bed
    """

rule compile_peaks:
    input:
        macs2="{outdir}/peaks/{sample}/macs2/{sample}_{genome}__peaks.broadPeak",
        epic2="{outdir}/peaks/{sample}/epic2/{sample}_{genome}__peaks_epic2.bed"
    output:
        peaks="{outdir}/peaks/{sample}/{sample}_{genome}__compiled_peaks.bed",
        rda="{outdir}/peaks/{sample}/{sample}_{genome}__compiled_peaks.rda"
    conda: helpers_dir + "/envs/compile_peaks.yaml"
    params:
        helpers_dir=helpers_dir,
        control=get_control,
        mode=get_mode,
        ouput_prefix="{outdir}/peaks/{sample}/{sample}_{genome}__compiled_peaks"
    log: "{outdir}/logs/compile_peaks/{sample}_{genome}__compile_peaks.log"
    shell: """
    (
    Rscript {params.helpers_dir}/scripts/compile_peaks.R {params.mode} {params.control} {input.macs2} {input.epic2} {params.ouput_prefix}
    ) &> {log}
    """

rule macs:
    input: unpack(input_test_callpeak)
    output: "{outdir}/peaks/{sample}/macs2/{sample}_{genome}__peaks.broadPeak"
    log: "{outdir}/logs/macs2/{sample}_{genome}__macs2.log"
    threads: 1
    conda: helpers_dir + "/envs/macs.yaml"
    params:
        prefix="{outdir}/peaks/{sample}/macs2/{sample}_{genome}_"
    shell: """
        (
        if [ {input.control} == {input.treatment} ]; then
            echo "No Control file detected -- running MACS2 without a control"
            macs2 callpeak --broad -t {input.treatment} -n {params.prefix}
        else
            echo "Control file detected -- running MACS2 with control"
            macs2 callpeak --broad -t {input.treatment} -c {input.control} -n {params.prefix}
        fi
        ) &> {log}
    """

rule epic:
    input:
        unpack(input_test_callpeak)
    output:
        "{outdir}/peaks/{sample}/epic2/{sample}_{genome}__peaks_epic2.bed"
    log: "{outdir}/logs/epic2/{sample}_{genome}__epic2_callpeaks.log"
    conda: helpers_dir + "/envs/epic.yaml"
    shell: """
        (
        if [ {input.control} == {input.treatment} ]; then
            echo "No Control file detected -- running EPIC2 without a control"
            epic2 -t {input.treatment} -gn {wildcards.genome} -o {output}
        else
            echo "Control file detected -- running EPIC2 with control"
            epic2 -t {input.treatment} -c {input.control} -gn {wildcards.genome} -o {output}
        fi
        ) &> {log}
    """

# This will be the rule that this executed if user provides bam file input
rule wrangle_bam:
    input: get_sample_bam
    output:
        bam = "{outdir}/wrangled_bam/{sample}/{sample}.{genome}.bam",
        bai = "{outdir}/wrangled_bam/{sample}/{sample}.{genome}.bam.bai"
    conda: helpers_dir + "/envs/bwa_mem.yaml"
    priority: 15
    log: "{outdir}/logs/wrangle_bam/{sample}_{genome}__wrangle_bam.log"
    params:
        bwa_extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        bwa_interleaved=pe_test_bwa,
        samblaster_extra=pe_test_samblaster,
        samtools_sort_extra="-O BAM"
    threads: 10
    shell: """
        (samtools view -h -@ {threads} -q 10 {input} | \
        samblaster {params.samblaster_extra}| \
        samtools sort {params.samtools_sort_extra} -@ {threads} -O bam -o {output.bam} - && \
        samtools index -@ {threads} {output.bam}) &> {log}
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
    input: check_type_fq
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
    input: find_sra_replicates
    output: temp("{outdir}/tmp/fastqs_merged/{sample}.{genome}__merged.fastq")
    log: "{outdir}/logs/merge_fastq/{sample}_{genome}__merge_fastq.log"
    shell: """
    (cat {input} > {output}) &> {log}
    """
    
    
rule cp_fq:
    input: find_fq
    output: temp("{outdir}/tmp/fastqs_cp/{sample}.{genome}__cp.fastq")
    log: "{outdir}/logs/cp_fastq/{sample}_{genome}__cp_fastq.log"
    shell: """
        (cp {input} {output}) &> {log}
    """
    

rule gunzip_fq:
    input: find_fq
    output: temp("{outdir}/tmp/fastqs_gunzip/{sample}.{genome}__gunzip.fastq")
    log: "{outdir}/logs/gunzip_fastq/{sample}_{genome}__gunzip_fastq.log"
    shell: """
        (cp {input} {output}.gz
        gunzip {output}.gz) &> {log}
    """


# If debugging, will use a less stable but much faster version of this rule which outputs a small fastq file
if debug:
    rule sra_to_fastq:
        # input: "{outdir}/tmp/sras/{sample}/{srr_acc}/{srr_acc}.sra"
        output: temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.fastq")
        conda: helpers_dir + "/envs/sratools.yaml"
        threads: 1
        log: "{outdir}/logs/sra_to_fastq/{sample}_{srr_acc}__sra_to_fastq_pe.log"
        params:
            output_directory="{outdir}/tmp/sras/{sample}/",
            fqdump="--skip-technical --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --split-3 "
        shell: """(
        fastq-dump -X 500000 {params.fqdump}-O {wildcards.outdir}/tmp/fastqs_raw/{wildcards.sample}/ {wildcards.srr_acc}
        cd {wildcards.outdir}/tmp/fastqs_raw/{wildcards.sample}/
        if test -f {wildcards.srr_acc}_2.fastq; then
            echo "Paired end -- interleaving"
            reformat.sh in1={wildcards.srr_acc}_1.fastq in2={wildcards.srr_acc}_2.fastq out={wildcards.srr_acc}.fastq overwrite=true
            rm {wildcards.srr_acc}_1.fastq && rm {wildcards.srr_acc}_2.fastq
        else
            echo "Single end -- finished!"
        fi
        ) &> {log}
        """
else:
    # This should always produced interleaved fq even if paired end. I don't like this solution, but it should hold.
    # reformat.sh (from BBTools) will interleave the paired-end files
    rule sra_to_fastq:
        # input: "{outdir}/tmp/sras/{sample}/{srr_acc}/{srr_acc}.sra"
        output: temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.fastq")
        conda: helpers_dir + "/envs/sratools.yaml"
        threads: 1
        log: "{outdir}/logs/sra_to_fastq/{sample}_{srr_acc}__sra_to_fastq_pe.log"
        params:
            output_directory="{outdir}/tmp/sras/{sample}/",
            fqdump="--skip-technical --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --split-3 "
        shell: """(
        cd {params.output_directory}
        fastq-dump {params.fqdump}-O ../../fastqs_raw/{wildcards.sample}/ {wildcards.srr_acc}
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
