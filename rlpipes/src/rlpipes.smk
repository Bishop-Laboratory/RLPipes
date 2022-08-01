########################################################################################################################
##############################################   Parse inputs    #######################################################
########################################################################################################################

import math
from os.path import expanduser


def rchop(s, suffix):
    """
    Python<3.9-friendly version of removesuffix()
    https://stackoverflow.com/questions/3663450/remove-substring-only-at-the-end-of-string
    """
    if suffix and s.endswith(suffix):
        return s[:-len(suffix)]
    return s

# Global Configs
src=config['src']
genome_home_dir = expanduser("~") + "/.rlpipes_genomes"
outdir=config['run_dir']
outdir= rchop(outdir, "/")
bwa_mem2=config['bwamem2']
macs3=config['macs3']
noreport=config['noreport']
useaws=config['useaws']

# Sample info
mode=config['mode']
paired_end=config['paired_end']
genome=config['genome']
sample_type=config['file_type']
sample_name=config['name']
eff_gen_size=config['eff_genome_size']

# Group info
groupby=None
if config['groupby'] is not None:
    groupby=config['groupby'].split(" ")
noexp=config['noexp']

# Refers to the experiment ID (e.g., basename of a fastq file, SRX)
# Used as primary name for file paths in this workflow.
sample=config['experiment']  
# Refers to raw seq data (e.g., fastq/bam file, SRR)
run=config['run']  
control=[ctr if ctr != "" else None for ctr in config['control']]

# Generate the output file names
report_html = expand("{outdir}/rlseq_report/{sample}_{genome}.html",zip,
            sample=sample, outdir=[outdir for i in range(len(sample))],genome=genome)
report_data = expand("{outdir}/rlseq_report/{sample}_{genome}.rda",zip,
            sample=sample, outdir=[outdir for i in range(len(sample))],genome=genome)
peaks_out = expand("{outdir}/peaks/{sample}_{genome}.broadPeak",zip,
            sample=sample, outdir=[outdir for i in range(len(sample))],genome=genome)
coverage_out = expand("{outdir}/coverage/{sample}_{genome}.bw",zip,
            sample=sample, outdir=[outdir for i in range(len(sample))],genome=genome)
bamstats_out = expand("{outdir}/bam_stats/{sample}_{genome}__bam_stats.txt",zip,
            sample=sample, outdir=[outdir for i in range(len(sample))],genome=genome)

# Get collater inputs
peaks=[outdir + '/peaks/' + elem + "_" + genome[idx] + ".broadPeak" for idx, elem in enumerate(sample) if mode[idx] not in ["RNA-Seq", "RNA-seq"]]
coverage=[outdir + '/coverage/' + elem + "_" + genome[idx] + ".bw" for idx, elem in enumerate(sample) if mode[idx] not in ["RNA-Seq", "RNA-seq"]]
bamstats=[outdir + '/bam_stats/' + elem + "_" + genome[idx] + "__bam_stats.txt" for idx, elem in enumerate(sample) if mode[idx] not in ["RNA-Seq", "RNA-seq"]]
quant=[outdir + '/quant/' + elem + "_" + genome[idx] + "/quant.sf" for idx, elem in enumerate(sample) if mode[idx] in ["RNA-Seq", "RNA-seq"]]
html=[outdir + '/rlseq_report/' + elem + "_" + genome[idx] + ".html" for idx, elem in enumerate(sample) if mode[idx] not in ["RNA-Seq", "RNA-seq"]]
rda=[outdir + '/rlseq_report/' + elem + "_" + genome[idx] + ".rda" for idx, elem in enumerate(sample) if mode[idx] not in ["RNA-Seq", "RNA-seq"]]
collate_inputs = html + rda + quant

# For testing the workflow using SRA
debug=config['debug']

# Select bwa type
# BWA MEM2 is still in development and has a particularly problematic habit of over-zealous RAM usage
# https://github.com/bwa-mem2/bwa-mem2/issues/118
# User needs to have the option to choose classic BWA
if bwa_mem2:
    bwa_cmd="bwa-mem2"
    bwa_ind_prefix="bwa_mem2_index/{genome}"
    bwa_location="bwa_mem2_index/{genome}.bwt.2bit.64"
    bwa_yaml = src + "/envs/bwamem2.yaml"
else:
    bwa_cmd="bwa"
    bwa_ind_prefix="bwa_index/{genome}"
    bwa_location="bwa_index/{genome}.sa"
    bwa_yaml = src + "/envs/bwa.yaml"

# Select MACS type
# MACS3 is still in development, but it is much faster
if macs3:
    macs_cmd="macs3"
    macs_yaml = src + "/envs/macs3.yaml"
else:
    macs_cmd="macs2"
    macs_yaml = src + "/envs/macs2.yaml"
    
    
# Set output
final_outputs = outdir + "/report.html"
# TODO: This is the part that needs to change when RLSeq is available
final_outputs = quant + peaks + coverage + bamstats
if not noreport:
    final_outputs = final_outputs + html + rda


########################################################################################################################
############################################   Helper Functions   ######################################################
########################################################################################################################
### Helper functions for pipeline ###
### For some reason, they have to be in the same file... ###
def find_sra_replicates(wildcards):
    srr_list = [run[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    replicates = expand("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.fastq",
           outdir=wildcards.outdir, srr_acc=srr_list, sample=wildcards.sample)
    return replicates

def find_fq(wildcards):
    fq = [run[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    return fq

def find_fq_pe(wildcards):
    fq1 = [re.sub('\\~.+', "", run[idx]) for idx, element in enumerate(sample) if element == wildcards.sample][0]
    fq2 = [re.sub('.+\\~', "", run[idx]) for idx, element in enumerate(sample) if element == wildcards.sample][0]
    return [
        fq1, fq2
    ]

def check_type_fq(wildcards):
    # class Wildcards:
    #   def __init__(self):
    #     self.sample=sample[0]
    # wildcards=Wildcards()
    file_type = [sample_type[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    if file_type == "fastq":
        fq = [run[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
        pe = [paired_end[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
        if pe:
            # CASE: paired-end fq
            return "{outdir}/tmp/fastqs_interleave/{sample}_{genome}__interleave.fastq"
            # TODO: Still need gz case...
        else:
            if fq[-3:] == ".gz":
                # CASE: Fastq GZ available
                return "{outdir}/tmp/fastqs_gunzip/{sample}_{genome}__gunzip.fastq"
            else:
                # CASE: normal Fastq
                return "{outdir}/tmp/fastqs_cp/{sample}_{genome}__cp.fastq"
    else:
        # CASE: Public data accession
        return "{outdir}/tmp/fastqs_merged/{sample}_{genome}__merged.fastq"
    

def pe_test_fastp(wildcards):
    pe = [paired_end[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    if pe:
        res="--interleaved_in "
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

def get_pe_bam(wildcards):
    pe = [paired_end[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    if pe:
        res="-f BAMPE "
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

def get_gensize(wildcards):
    return [eff_gen_size[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    
def isdebug(wildcards):
    if debug:
        param="-X 500000 "
    else:
        param=""
    return param


def choose_bam_type(wildcards):
    """
    Switch to the "wrangled_bam" directory if user supplied a bam file.
    """
    st_now = [sample_type[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    bam_type = "/bam/"
    if st_now == "bam":
        bam_type = "/wrangled_bam/"
    return [
        wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.genome + ".bam",
        wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.genome + ".bam.bai"
    ]


def get_report_inputs(wildcards):
    mode_now = [mode[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    st_now = [sample_type[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    return_dict = {
        'peaks': wildcards.outdir + '/peaks/' + wildcards.sample + "_" + wildcards.genome + ".broadPeak",
        'coverage': wildcards.outdir + '/coverage/' + wildcards.sample + "_" + wildcards.genome + ".bw",
        'bam_stats': wildcards.outdir + '/bam_stats/' + wildcards.sample + "_" + wildcards.genome + "__bam_stats.txt",
    }
    return return_dict

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
            inpt = [sample[idx] for idx, element in enumerate(sample) if element == inpt][0]
    return_dict = {
            'treatment': wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.genome + ".bam",
            'control': wildcards.outdir + bam_type + inpt + "/" + inpt + "_" + wildcards.genome + ".bam",
            'index_treatment': wildcards.outdir + bam_type + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.genome + ".bam.bai",
            'index_control': wildcards.outdir + bam_type + inpt + "/" + inpt + "_" + wildcards.genome + ".bam.bai"
        }
    return return_dict

def test_pe(wildcards):
  return [paired_end[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]

def get_fq_salmon(wildcards):
  pe = test_pe(wildcards)
  if pe:
      res="-1 " + wildcards.outdir + "/tmp/fastqs_prepped/" + wildcards.sample + "_" + wildcards.genome +\
      ".R1.fastq -2 " + wildcards.outdir + "/tmp/fastqs_prepped/" + wildcards.sample + "_" +\
      wildcards.genome + ".R2.fastq"
  else:
      res="-r " + wildcards.outdir + "/tmp/fastqs_prepped/" + wildcards.sample + "_" + wildcards.genome + ".R1.fastq"
  return res

########################################################################################################################
##############################################   Main pipeline    ######################################################
########################################################################################################################

rule output:
    input: final_outputs
    
    
rule rlseq_collate:
    input: collate_inputs
    output:
        src="{outdir}/report.html"
    shell:
        """
        touch {output.src}
        """
        
rule salmon_quant:
  input: 
    fq="{outdir}/tmp/fastqs_prepped/{sample}_{genome}.R1.fastq",
    ind=genome_home_dir + "/{genome}/salmon_index/versionInfo.json"
  output: "{outdir}/quant/{sample}_{genome}/quant.sf"
  params:
    pe_fq=get_fq_salmon,
    ind= genome_home_dir + "/{genome}/salmon_index",
    outdir="{outdir}/quant/{sample}_{genome}"
  log: "{outdir}/logs/quant/{sample}_{genome}__quant.log"
  conda: src + "/envs/salmon.yaml"
  threads: 40
  priority: 11
  shell:
    """
    (
      salmon quant -i {params.ind} -l A {params.pe_fq} --validateMappings -o {params.outdir} -p {threads}
    ) &> {log}
    
    """

rule cleanup_fq:
  input: "{outdir}/tmp/fastqs_trimmed/{sample}_{genome}__trimmed.fastq"
  output: temp("{outdir}/tmp/fastqs_prepped/{sample}_{genome}.R1.fastq")
  conda: src + "/envs/sratools.yaml"
  log: "{outdir}/logs/cleanup_fq/{sample}.{genome}__cleanup_fq.log"
  params:
    ispe=test_pe,
    outdir="{outdir}/tmp/fastqs_prepped/"
  shell: """
  (
    if [ {params.ispe} == "True" ]; then
      reformat.sh ow=t in={input} out1={params.outdir}/{wildcards.sample}_{wildcards.genome}.R1.fastq \
      out2={params.outdir}/{wildcards.sample}_{wildcards.genome}.R2.fastq
    else
      mv {input} {output}
    fi
  ) &> {log}
  """
  

rule salmon_index:
    input:
        fa=genome_home_dir + "/{genome}/{genome}.ensTx.fa"
    output:
        vi=genome_home_dir + "/{genome}/salmon_index/versionInfo.json"
    params:
        prefix=genome_home_dir + "/{genome}/salmon_index"
    threads: 50
    conda: src + "/envs/salmon.yaml"
    log: genome_home_dir + "/{genome}/salmon_index/{genome}__salmon_index.log"
    shell: """
    (
        salmon index -t {input} -i {params.prefix} -p {threads}
    ) &> {log}
        
    """

rule prep_txome_fa:
    input:
        fa=genome_home_dir + "/{genome}/{genome}.fa",
        gtf=genome_home_dir + "/{genome}/{genome}.ensGene.gtf"
    output:
        genome_home_dir + "/{genome}/{genome}.ensTx.fa"
    params:
        pref=genome_home_dir + "/{genome}/rsem/{genome}",
        rsemdir=genome_home_dir + "/{genome}/rsem"
    conda: src + "/envs/rsem.yaml"
    log: genome_home_dir + "/logs/prep_txome/{genome}__prep_txome_fa.log"
    shell: """
    (
        if [ -d {params.rsemdir} ]; then rm -Rf {params.rsemdir}; fi
        mkdir {params.rsemdir}
        rsem-prepare-reference --gtf {input.gtf} {input.fa} {params.pref}
        mv {params.pref}.transcripts.fa {output}
        rm -rf {params.rsemdir}
    ) &> {log}
        
    """

rule download_gtf:
    output:
        genome_home_dir + "/{genome}/{genome}.ensGene.gtf"
    log: genome_home_dir + "/logs/download_gtf/{genome}__download_gtf.log"
    shell: """
        wget -O {output}.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/genes/{wildcards.genome}.ensGene.gtf.gz
        gunzip {output}.gz &> {log}
    """

rule rlseq:
    input: unpack(get_report_inputs)
    output:
        html="{outdir}/rlseq_report/{sample}_{genome}.html",
        data="{outdir}/rlseq_report/{sample}_{genome}.rda"
    params:
        src=src,
        configs = "{outdir}/config.json",
    log: "{outdir}/logs/rlseq_report/{sample}_{genome}__rlseq.log"
    conda: src + "/envs/rlseq.yaml"
    script: src + "/scripts/RLSeq.R"
    
rule calculate_coverage:
    input: choose_bam_type
    output: "{outdir}/coverage/{sample}_{genome}.bw"
    threads: 10
    conda: src + "/envs/deeptools.yaml"
    log: "{outdir}/logs/coverage/{sample}_{genome}__coverage.log"
    params:
        extra="--minMappingQuality 20 --ignoreDuplicates"
    shell: """
        (bamCoverage -b {input[0]} -p {threads} {params.extra} -o {output}) &> {log}
    """

rule bam_stats:
    input: choose_bam_type
    output: "{outdir}/bam_stats/{sample}_{genome}__bam_stats.txt"
    threads: 4
    conda: src + "/envs/samtools.yaml"
    log: "{outdir}/logs/bam_stats/{sample}_{genome}__bam_stats.log"
    shell: "(samtools flagstat -@ {threads} {input[0]} > {output}) &> {log}"


rule macs_callpeak:
    input: unpack(input_test_callpeak)
    output: "{outdir}/peaks/{sample}_{genome}.broadPeak"
    log: "{outdir}/logs/peaks/{sample}_{genome}__macs.log"
    threads: 1
    conda: macs_yaml
    params:
        prefix="{outdir}/peaks/{sample}_{genome}_",
        macsout="{outdir}/peaks/{sample}_{genome}__peaks.broadPeak",
        gensize=get_gensize,
        macs_cmd=macs_cmd,
        pe_bam=get_pe_bam
    shell: """
        (
        if [ {input.control} == {input.treatment} ]; then
            echo {params.macs_cmd}
            echo "No Control file detected -- running MACS without a control"
            {params.macs_cmd} callpeak --broad {params.pe_bam}-g {params.gensize} \
            -t {input.treatment} -n {params.prefix} || true
        else
            echo "Control file detected -- running MACS with control"
            echo {params.macs_cmd}
            {params.macs_cmd} callpeak --broad {params.pe_bam}-g {params.gensize} \
            -t {input.treatment} -c {input.control} -n {params.prefix} || true
        fi
        mv {params.macsout} {output} || (touch {output} && touch {output}.failed)
        ) &> {log}
    """


# This will be the rule that this executed if user provides bam file input
rule wrangle_bam:
    input: get_sample_bam
    output:
        bam = "{outdir}/wrangled_bam/{sample}/{sample}_{genome}.bam",
        bai = "{outdir}/wrangled_bam/{sample}/{sample}_{genome}.bam.bai"
    conda: src + "/envs/samtools.yaml"
    priority: 15
    log: "{outdir}/logs/wrangle_bam/{sample}_{genome}__wrangle_bam.log"
    params:
        bwa_extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        bwa_interleaved=pe_test_bwa,
        samtools_sort_extra="-O BAM"
    threads: 30
    shell: """
        (   
            samtools view -h -@ {threads} -q 10 {input} | \
            samtools sort {params.samtools_sort_extra} -@ {threads} -O bam -o {output.bam} - && \
            samtools index -@ {threads} {output.bam}
        ) &> {log}
     """


rule index_bam:
    input: 
        bam="{outdir}/bam/{sample}/{sample}_{genome}.bam"
    output: 
        bai="{outdir}/bam/{sample}/{sample}_{genome}.bam.bai"
    conda: src + "/envs/samtools.yaml"
    threads: 10
    log: "{outdir}/logs/index_bam/{sample}_{genome}__index_bam.log"
    shell: """
        samtools index -@ {threads} {input.bam}
    """


rule bwa_mem:
    input:
        bwa_index_done=ancient(genome_home_dir + "/{genome}/" + bwa_location),
        reads=ancient("{outdir}/tmp/fastqs_trimmed/{sample}_{genome}__trimmed.fastq")
    output:
        bam="{outdir}/bam/{sample}/{sample}_{genome}.bam"
    conda: bwa_yaml
    priority: 15
    log: "{outdir}/logs/bwa/{sample}_{genome}__bwa_mem.log"
    params:
        index=genome_home_dir + "/{genome}/" + bwa_ind_prefix,
        bwa_extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        bwa_interleaved=pe_test_bwa,
        samtools_sort_extra="-O BAM",
        bwa_cmd=bwa_cmd
    threads: 40
    shell: """
        ({params.bwa_cmd} mem -t {threads} {params.bwa_extra} {params.bwa_interleaved}{params.index} {input.reads} | \
        samtools view -q 10 -b -@ {threads} - | \
        samtools sort {params.samtools_sort_extra} -@ {threads} -o {output.bam} -) &> {log}
     """


rule bwa_index:
    input:
        genome_home_dir + "/{genome}/{genome}.fa"
    output:
        genome_home_dir + "/{genome}/" + bwa_location,
    params:
        prefix=genome_home_dir + "/{genome}/" + bwa_ind_prefix,
        bwa_cmd=bwa_cmd
    conda: bwa_yaml
    log:
        genome_home_dir + "/{genome}/" + bwa_ind_prefix + "_bwa_index.log"
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
        trimmed=temp("{outdir}/tmp/fastqs_trimmed/{sample}_{genome}__trimmed.fastq"),
        json="{outdir}/fastq_stats/{sample}_{genome}__fastq_stats.json"
    conda: src + "/envs/fastp.yaml"
    log: "{outdir}/logs/fastp/{sample}_{genome}__fastp_pe.log"
    priority: 10
    params:
        extra=pe_test_fastp
    threads: 4
    shell: """
    (fastp -i {input} --stdout {params.extra}-w {threads} -h /dev/null -j {output.json} > {output} ) &> {log}
    """


rule cp_fq:
    input: find_fq
    output: temp("{outdir}/tmp/fastqs_cp/{sample}_{genome}__cp.fastq")
    log: "{outdir}/logs/cp_fastq/{sample}_{genome}__cp_fastq.log"
    shell: """
        (cp {input} {output}) &> {log}
    """
    

rule gunzip_fq:
    input: find_fq
    output: temp("{outdir}/tmp/fastqs_gunzip/{sample}_{genome}__gunzip.fastq")
    log: "{outdir}/logs/gunzip_fastq/{sample}_{genome}__gunzip_fastq.log"
    shell: """
        (cp {input} {output}.gz
        gunzip {output}.gz) &> {log}
    """
    

rule fq_interleave:
    input: 
        r1="{outdir}/tmp/fastqs_repaired/{sample}_{genome}__repair.R1.fastq",
        r2="{outdir}/tmp/fastqs_repaired/{sample}_{genome}__repair.R2.fastq"
    output: temp("{outdir}/tmp/fastqs_interleave/{sample}_{genome}__interleave.fastq")
    log: "{outdir}/logs/fastqs_interleave/{sample}_{genome}__interleave_fastq.log"
    conda: src + "/envs/sratools.yaml"
    shell:"""
    (
        echo "Paired end -- interleaving"
        reformat.sh in1={input.r1} in2={input.r2} out={output} overwrite=true
    ) &> {log}
    """


rule fq_repair:
    input: find_fq_pe
    output:
        r1=temp("{outdir}/tmp/fastqs_repaired/{sample}_{genome}__repair.R1.fastq"),
        r2=temp("{outdir}/tmp/fastqs_repaired/{sample}_{genome}__repair.R2.fastq")
    log: "{outdir}/logs/fastqs_repair/{sample}_{genome}__repair_fastq.log"
    conda: src + "/envs/sratools.yaml"
    shell: """
    (
        echo "Repairing mates"
        repair.sh in1={input[0]} in2={input[1]} out1={output.r1} out2={output.r2} outs=/dev/null repair
    ) &> {log}
    """


rule merge_replicate_reads:
    input: find_sra_replicates
    output: temp("{outdir}/tmp/fastqs_merged/{sample}_{genome}__merged.fastq")
    log: "{outdir}/logs/merge_fastq/{sample}_{genome}__merge_fastq.log"
    shell: """
    (cat {input} > {output}) &> {log}
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
if not config['useaws']:
  rule download_sra_prefetch:
      output: temp("{outdir}/tmp/sras/{sample}/{srr_acc}/{srr_acc}.sra")
      conda: src + "/envs/sratools.yaml"
      log: "{outdir}/logs/download_sra_prefetch/{sample}__{srr_acc}__download_sra_prefetch.log"
      params:
          output_directory = "{outdir}/tmp/sras/{sample}/"
      shell: """
              (
              cd {params.output_directory}
              prefetch {wildcards.srr_acc} -f yes
              ) &> {log}
              """
else:
  rule download_sra_aws:
      output: temp("{outdir}/tmp/sras/{sample}/{srr_acc}/{srr_acc}.sra")
      log: "{outdir}/logs/download_sra_aws/{sample}__{srr_acc}__download_sra_aws.log"
      conda: src + "/envs/awscli.yaml"
      shell: """
              (
              echo "Downloading from AWS..."
              aws s3 cp s3://sra-pub-run-odp/sra/{wildcards.srr_acc}/{wildcards.srr_acc} {output} --no-sign-request --no-progress
              ) &> {log}
              """
