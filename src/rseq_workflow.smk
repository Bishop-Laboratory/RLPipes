########################################################################################################################
##############################################   Parse inputs    #######################################################
########################################################################################################################

import math
from os.path import expanduser

# Global Configs
src=config['src'][0]
genome_home_dir = expanduser("~") + "/.rseq_genomes"
outdir=config['outdir'][0]
outdir=outdir.removesuffix("/")
bwa_mem2=config['bwamem2'][0]
macs3=config['macs3'][0]

# Sample info
mode=config['mode']
paired_end=config['paired_end']
genome=config['genome']
sample_type=config['file_type']
sample_name=config['name']
# Refers to the experiment ID (e.g., basename of a fastq file, SRX)
# Used as primary name for file paths in this workflow.
sample=config['experiment']  
# Refers to raw seq data (e.g., fastq/bam file, SRR)
run=config['run']  
control=[ctr if ctr != "" else None for ctr in config['control']]

# Generate the output file names
report_html = expand("{outdir}/RSeq_report/{sample}_{genome}.html",zip,
            sample=sample, outdir=[outdir for i in range(len(sample))],genome=genome)
report_data = expand("{outdir}/RSeq_report/{sample}_{genome}.rda",zip,
            sample=sample, outdir=[outdir for i in range(len(sample))],genome=genome)

# For testing the workflow using SRA
debug=False

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
    
# Select MACS type
# MACS3 is still in development, but it is much faster
if macs3:
    macs_cmd="macs3"
    macs_yaml = src + "/envs/macs3.yaml"
else:
    macs_cmd="macs2"
    macs_yaml = src + "/envs/macs2.yaml"

### Helper functions for pipeline ###
def find_sra_replicates(wildcards):
    srr_list = [run[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    replicates = expand("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.fastq",
           outdir=wildcards.outdir, srr_acc=srr_list, sample=wildcards.sample)
    return replicates

def find_fq(wildcards):
    fq = [run[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    return fq

def check_type_fq(wildcards):
    file_type = [sample_type[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    if file_type == "fastq":
        fq = [run[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0][0]
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
    return [control[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]

def get_sample_bam(wildcards):
    return [run[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]

def isdebug(wildcards):
    if debug:
        param="-X 500000 "
    else:
        param=""
    return param

def input_test_callpeak(wildcards):
    inpt = [control[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    st_now = [sample_type[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    # This is a hack that changes the expected bam location to switch between the bwamem rule and bam wrangle rule
    if st_now != "bam":
        bam_type = "/bam/"
    else:
        bam_type = "/wrangled_bam/"

    # Set the expected samples for peak calling
    if inpt is None:
        inpt = wildcards.sample
    else:
        if st_now != "public":
            # If not public, input will be a file. Need input sample name instead.
            inpt = [sample[idx] for idx, element in enumerate(samples) if element[0] == inpt][0]
    return_dict = {
            'treatment': wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam",
            'control': wildcards.outdir + bam_type + inpt + "/" + inpt + "." + wildcards.genome + ".bam",
            'index_treatment': wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam.bai",
            'index_control': wildcards.outdir + bam_type + inpt + "/" + inpt + "." + wildcards.genome + ".bam.bai"
        }
    return return_dict

def get_report_inputs(wildcards):
    return_dict = {
        'peaks': wildcards.outdir + '/peaks/' + wildcards.sample + "_" + wildcards.genome + ".broadPeak",
        'coverage': wildcards.outdir + '/coverage/' + wildcards.sample + "_" + wildcards.genome + ".bw",
        'bam_stats': wildcards.outdir + '/bam_stats/' + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.genome + "__bam_stats.txt"
    }
    st_now = [sample_type[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    if st_now in ['fastq', 'public']:
        return_dict['read_qc_data'] = wildcards.outdir + '/QC/fastq/json/' + wildcards.sample + "." + wildcards.genome + ".json"
    return return_dict

def choose_bam_type(wildcards):
    """
    Switch to the "wrangled_bam" directory if user supplied a bam file.
    """
    st_now = [sample_type[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    bam_type = "/bam/"
    if st_now == "bam":
        bam_type = "/wrangled_bam/"
    return wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "." + wildcards.genome + ".bam"


#######################################################################################################################
##############################################   Main pipeline    ######################################################
########################################################################################################################

rule output:
    input:
        html=report_html,
        data=report_data
        
        
rule RSeqR:
    input: unpack(get_report_inputs)
    output:
        html="{outdir}/RSeq_report/{sample}_{genome}.html",
        data="{outdir}/RSeq_report/{sample}_{genome}.rda"
    # conda: src + "/envs/rseqr.yaml"  TODO: Fix for production once RSeqR is public
    params:
        src=src,
        configs = "{outdir}/config.json",
    log: "{outdir}/logs/prepare_report/{sample}_{genome}__prepare_report.log"
    script: src + "scripts/runRSeqR.R"
    
    
rule calculate_coverage:
    input: choose_bam_type
    output: "{outdir}/coverage/{sample}_{genome}.bw"
    conda: src + "/envs/deeptools.yaml"
    log: "{outdir}/logs/coverage/{sample}_{genome}__coverage.log"
    threads: 5
    params:
        extra="--minMappingQuality 20 --binSize 10"
    shell: """
        (bamCoverage -b {input} -p {threads} {params.extra} -o {output}) &> {log}
    """

rule bam_stats:
    input: choose_bam_type
    output: "{outdir}/bam_stats/{sample}/{sample}_{genome}__bam_stats.txt"
    threads: 4
    conda: src + "/envs/samtools.yaml"
    log: "{outdir}/logs/bam_stats/{sample}_{genome}__bam_stats.log"
    shell: "(samtools flagstat -@ {threads} {input} > {output}) &> {log}"


rule macs_callpeak:
    input: unpack(input_test_callpeak)
    output: "{outdir}/peaks/{sample}_{genome}.broadPeak"
    log: "{outdir}/logs/peaks/{sample}_{genome}__macs2.log"
    conda: macs_yaml
    params:
        prefix="{outdir}/peaks/{sample}/{sample}_{genome}_",
        macsout="{outdir}/peaks/{sample}_{genome}__peaks.broadPeak",
        macs_cmd=macs_cmd
    shell: """
        (
        if [ {input.control} == {input.treatment} ]; then
            echo "No Control file detected -- running MACS2 without a control"
            {params.macs_cmd} callpeak --broad -t {input.treatment} -n {params.prefix}
        else
            echo "Control file detected -- running MACS2 with control"
            {params.macs_cmd} callpeak --broad -t {input.treatment} -c {input.control} -n {params.prefix}
        fi
        mv {params.macsout} {output}
        ) &> {log}
    """


# This will be the rule that this executed if user provides bam file input
rule wrangle_bam:
    input: get_sample_bam
    output:
        bam = "{outdir}/wrangled_bam/{sample}/{sample}.{genome}.bam",
        bai = "{outdir}/wrangled_bam/{sample}/{sample}.{genome}.bam.bai"
    conda: src + "/envs/bwa_mem.yaml"
    threads: 30
    log: "{outdir}/logs/wrangle_bam/{sample}_{genome}__wrangle_bam.log"
    params:
        bwa_extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        bwa_interleaved=pe_test_bwa,
        samblaster_extra=pe_test_samblaster,
        samtools_sort_extra="-O BAM"
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
            genome_home_dir + "/{genome}/bwa_index/{genome}.amb"
        ],
        reads="{outdir}/tmp/fastqs_trimmed/{sample}.{genome}__trimmed.fastq"
    output:
        bam="{outdir}/bam/{sample}/{sample}.{genome}.bam",
        bai="{outdir}/bam/{sample}/{sample}.{genome}.bam.bai"
    conda: src + "/envs/bwa_mem.yaml"
    threads: 30
    log: "{outdir}/logs/bwa/{sample}_{genome}__bwa_mem.log"
    params:
        index=genome_home_dir + "/{genome}/bwa_index/{genome}",
        bwa_extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        bwa_interleaved=pe_test_bwa,
        samblaster_extra=pe_test_samblaster,
        samtools_sort_extra="-O BAM",
        bwa_cmd=bwa_cmd
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
    conda: src + "/envs/bwa.yaml"
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
    conda: src + "/envs/fastp.yaml"
    log: "{outdir}/logs/fastp/{sample}.{genome}__fastp_pe.log"
    params:
        extra=pe_test_fastp
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


# This should always produced interleaved fq even if paired end. I don't like this solution, but it should hold.
# reformat.sh (from BBTools) will interleave the paired-end files
rule sra_to_fastq:
    input: "{outdir}/tmp/sras/{sample}/{srr_acc}/{srr_acc}.sra"
    output: temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.fastq")
    conda: src + "/envs/sratools.yaml"
    log: "{outdir}/logs/sra_to_fastq/{sample}_{srr_acc}__sra_to_fastq_pe.log"
    params:
        output_directory="{outdir}/tmp/sras/{sample}/",
        fqdump="--skip-technical --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --split-3 ",
        debug=isdebug
    shell: """(
    cd {params.output_directory}
    fastq-dump {params.debug}{params.fqdump}-O ../../fastqs_raw/{wildcards.sample}/ {wildcards.srr_acc}
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
    conda: src + "/envs/sratools.yaml"
    log: "{outdir}/logs/download_sra/{sample}__{srr_acc}__download_sra.log"
    params:
        output_directory = "{outdir}/tmp/sras/{sample}/"
    shell: """
            (
            cd {params.output_directory}
            prefetch {wildcards.srr_acc} -f yes
            ) &> {log}
            """
