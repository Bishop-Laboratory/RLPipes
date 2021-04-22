########################################################################################################################
##############################################   Parse inputs    #######################################################
########################################################################################################################

# Global Configs
helpers_dir=config['helpers_dir']
genome_home_dir=config['genome_home_dir'][0]
cores=config['threads'][0]
outdir=config['outdir'][0]
outdir=outdir.strip("/")

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
controls=[ctr.split(",") for ctr in config['control'] if ctr != "None"]

bam_index_output = expand("{outdir}/bam/{sample}/{sample}.{genome}.bam.bai", zip,
                        sample=sample, outdir=[outdir for i in range(len(sample))], genome=genome)
coverage_output = expand("{outdir}/coverage/{sample}/{sample}.{genome}.bw", zip,
                        sample=sample, outdir=[outdir for i in range(len(sample))], genome=genome)



#######################################################################################################################
##############################################   Main pipeline    ######################################################
########################################################################################################################
def find_replicates(wildcards):
    srr_list = [experiments[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
    replicates = expand("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.fastq",
           outdir=wildcards.outdir, srr_acc=srr_list, sample=wildcards.sample)
    return replicates

def pe_test_fastp(wildcards):
    pe = [paired_end[idx] for idx, element in enumerate(sample) if element == wildcards.sample]
    if pe[0]:
        res="--interleaved_in "
    else:
        res=""
    return res

def pe_test_samblaster(wildcards):
    pe = [paired_end[idx] for idx, element in enumerate(sample) if element == wildcards.sample]
    if not pe[0]:
        res="--ignoreUnmated "
    else:
        res=""
    return res

def pe_test_bwa(wildcards):
    pe = [paired_end[idx] for idx, element in enumerate(sample) if element == wildcards.sample]
    if pe[0]:
        res="-p "
    else:
        res=""
    return res

rule output:
    input:
        bam_index_output,
        # coverage_output

rule download_fasta:
    output:
        genome_home_dir + "/{genome}/{genome}.fa"
    params:
          prefix=genome_home_dir + "/{genome}/bwa_index/{genome}",
          check=outdir + "/logs/download_fasta/{genome}__download_fasta.log"
    shell: """
        (mkdir -p {params.prefix}
        wget -O {output}.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz
        gunzip {output}.gz && echo "Indexing BAM at {params.prefix} ... This may take some time ...") &> {params.check}
    """

rule bwa2_index:
    input:
        genome_home_dir + "/{genome}/{genome}.fa"
    output:
        genome_home_dir + "/{genome}/bwa_index/{genome}.ann",
        genome_home_dir + "/{genome}/bwa_index/{genome}.bwt.2bit.64",
        genome_home_dir + "/{genome}/bwa_index/{genome}.pac",
        genome_home_dir + "/{genome}/bwa_index/{genome}.amb"
    params:
        prefix=genome_home_dir + "/{genome}/bwa_index/{genome}"
    conda: helpers_dir + "/envs/bwa.yaml"
    log:
        genome_home_dir + "/{genome}/bwa_index/{genome}_bwa_index.log"
    shell:"""
        (bwa-mem2 index -p {params.prefix} {input}) &> {log}
    """

# TODO: Figure out how to use the pipes
# TODO: Probably need to specify this version of prefetch and/or find alternative to it...
# TODO: Retry if fails due to network error
rule download_sra:
    output: temp("{outdir}/tmp/sras/{sample}/{srr_acc}/{srr_acc}.sra")
    conda: helpers_dir + "/envs/sratools.yaml"
    log: "{outdir}/logs/download_sra/{sample}__{srr_acc}__download_sra.log"
    params:
        output_directory="{outdir}/tmp/sras/{sample}/"
    threads: cores*.2 - 1
    shell: """
    (
    cd {params.output_directory}
    prefetch {wildcards.srr_acc} -f yes
    ) &> {log}
    """

# rule sra_to_fastq:
#     input: "{outdir}/tmp/sras/{sample}/{srr_acc}/{srr_acc}.sra"
#     output: temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.fastq")
#     conda: helpers_dir + "/envs/parallel_fq_dump.yaml"
#     threads: 8
#     log: "{outdir}/logs/sra_to_fastq/{sample}_{srr_acc}__sra_to_fastq_pe.log"
#     params:
#         output_directory="{outdir}/tmp/sras/{sample}/"
#     shell: """(
#     cd {params.output_directory}
#     parallel-fastq-dump -t {threads} --split-3 -Z -s {wildcards.srr_acc} --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' > ../../fastqs_raw/{wildcards.sample}/{wildcards.srr_acc}.fastq
#     ) &> {log}
#     """

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

# # TODO: Really should be a better way than this...
# # TODO: Still need to find a situation with un-mated reads -- may need a custom solution for that
# # Okay, so here's the problem: fastq-dump will create
# rule sra_to_fastq:
#     input: "{outdir}/tmp/sras/{sample}/{srr_acc}/{srr_acc}.sra"
#     output: temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.fastq")
#     conda: helpers_dir + "/envs/sratools.yaml"
#     threads: 1
#     log: "{outdir}/logs/sra_to_fastq/{sample}_{srr_acc}__sra_to_fastq_pe.log"
#     params:
#         output_directory="{outdir}/tmp/sras/{sample}/"
#     shell: """(
#     cd {params.output_directory}
#     fast-dump -X 100000 --split-3 -O ../../fastqs_raw/{wildcards.sample}/ --skip-technical --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' {wildcards.srr_acc}
#     ) &> {log}
#     """

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

# TODO: Try BWA MEM2 (Has stdin option and 1.5-3X speed increase)
rule bwa_mem2:
    input:
        bwa_index_done=[
              genome_home_dir + "/{genome}/bwa_index/{genome}.ann",
              genome_home_dir + "/{genome}/bwa_index/{genome}.bwt.2bit.64",
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
        samtools_sort_extra="-O BAM"
    threads: 10
    shell: """
        (bwa-mem2 mem -t {threads} {params.bwa_extra} {params.bwa_interleaved}{params.index} {input.reads} | \
        samblaster {params.samblaster_extra}| \
        samtools view -q 10 -b -@ {threads} - | \
        samtools sort {params.samtools_sort_extra} -@ {threads} -o {output.bam} - && \
        samtools index -@ {threads} {output.bam}) &> {log}
     """







#
#
#
# if sample_type != "bam" and sample_type != "bigWig" and sample_type != "bedGraph":
#     if sample_type == "public":
#
#
#     if paired_end:
#         if sample_type == "public":
#             # TODO: fastq merging only works with public samples at the moment. This should also work with local samples.
#             rule sra_to_fastq_pe:
#                 input:
#                     "{outdir}/tmp/sras/{sample}/{srr_acc}.sra"
#                 output:
#                     temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}_1.fastq"),
#                     temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}_2.fastq")
#                 threads: threads
#                 log: "{outdir}/logs/{sample}_{srr_acc}__sra_to_fastq_pe.log"
#                 shell:
#                     "(parallel-fastq-dump -t {threads} --split-files -O" \
#                     + " {wildcards.outdir}/tmp/fastqs_raw/{wildcards.sample}/ -s {input}) &> {log}"
#
#         rule merge_fastq_pe_experiment:
#             input:
#                 R1=merge_input_experiment_1,
#                 R2=merge_input_experiment_2
#             output:
#                 R1=temp("{outdir}/tmp/fastqs_dup/{sample}_experiment_R1.fastq"),
#                 R2=temp("{outdir}/tmp/fastqs_dup/{sample}_experiment_R2.fastq")
#             log: "{outdir}/logs/{sample}_merge_fastq_pe_experiment.log"
#             shell: """
#                 ( cat {input.R1} > {output.R1}
#                 cat {input.R2} > {output.R2} ) &> {log}
#             """
#
#         if controls != "None":
#             rule merge_fastq_pe_control:
#                 input:
#                     R1=merge_input_control_1,
#                     R2=merge_input_control_2
#                 output:
#                     R1=temp("{outdir}/tmp/fastqs_dup/{sample}_control_R1.fastq"),
#                     R2=temp("{outdir}/tmp/fastqs_dup/{sample}_control_R2.fastq")
#                 log: "{outdir}/logs/{sample}_merge_fastq_pe_control.log"
#                 shell: """
#                     ( cat {input.R1} > {output.R1}
#                     cat {input.R2} > {output.R2} ) &> {log}
#                 """
#
#         if not no_fastp:
#             rule fastp_pe:
#                 input:
#                     sample=["{outdir}/tmp/fastqs_dup/{sample}_{exp_type}_R1.fastq",
#                             "{outdir}/tmp/fastqs_dup/{sample}_{exp_type}_R2.fastq"]
#                 output:
#                     trimmed=[temp("{outdir}/fastqs/{sample}_{exp_type}_R1.fastq"),
#                              temp("{outdir}/fastqs/{sample}_{exp_type}_R2.fastq")],
#                     html="{outdir}/QC/fastq/{sample}.{exp_type}.html",
#                     json="{outdir}/QC/fastq/{sample}.{exp_type}.json"
#                 log: "{outdir}/logs/{sample}_{exp_type}__fastp_pe.log"
#                 params:
#                     extra=""
#                 threads: threads
#                 wrapper:
#                     "0.63.0/bio/fastp"
#
#     else:
#         if sample_type == "public":
#             rule sra_to_fastq_se:
#                 input: "{outdir}/tmp/sras/{sample}/{srr_acc}.sra"
#                 output: temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}_1.fastq")
#                 threads: threads
#                 log: "{outdir}/logs/{sample}_{srr_acc}__sra_to_fastq_se.log"
#                 shell:
#                     "(parallel-fastq-dump -t {threads} --split-files -O" \
#                     + " {wildcards.outdir}/tmp/fastqs_raw/{wildcards.sample}/ -s {input}) &> {log}"
#
#         rule merge_fastq_se_experiment:
#             input: merge_input_experiment
#             output: temp("{outdir}/tmp/fastqs_dup/{sample}_experiment.fastq")
#             log: "{outdir}/logs/{sample}_merge_fastq_se_experiment.log"
#             shell: "(cat {input} > {output}) &> {log}"
#
#         if controls != "None":
#             rule merge_fastq_se_control:
#                 input: merge_input_control
#                 output: temp("{outdir}/tmp/fastqs_dup/{sample}_control.fastq")
#                 log: "{outdir}/logs/{sample}_merge_fastq_se_control.log"
#                 shell: "(cat {input} > {output}) &> {log}"
#
#         if not no_fastp:
#             rule fastp_se:
#                 input:
#                      sample=["{outdir}/tmp/fastqs_dup/{sample}_{exp_type}.fastq"]
#                 output:
#                     trimmed=temp("{outdir}/fastqs/{sample}_{exp_type}.fastq"),
#                     html="{outdir}/QC/fastq/{sample}.{exp_type}.html",
#                     json="{outdir}/QC/fastq/{sample}.{exp_type}.json"
#                 log: "{outdir}/logs/{sample}_{exp_type}__fastp_se.log"
#                 params:
#                     extra=""
#                 threads: threads
#                 wrapper:
#                     "0.63.0/bio/fastp"
#
#     rule download_fasta:
#         output:
#             genome_home_dir + "/{genome}/{genome}.fa"
#         params:
#               prefix=genome_home_dir + "/{genome}/bwa_index/{genome}",
#               check=outdir + "/logs/" + sample_name + "_{genome}__download_fasta.log"
#         shell: """
#             (mkdir -p {params.prefix}
#             wget -O {output}.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz
#             gunzip {output}.gz && echo "Indexing BAM at {params.prefix} ... This may take some time ...") &> {params.check}
#         """
#
#     rule download_homer_anno:
#         output:
#             gene_anno=genome_home_dir + "/{genome}/homer_anno.txt"
#         params:
#             out_dir=genome_home_dir + "/{genome}/homer_anno",
#             check=outdir + "/logs/" + sample_name + "_{genome}__download_homer_anno.log"
#         shell: """
#             (wget -O {params.out_dir}.zip http://homer.ucsd.edu/homer/data/genomes/{wildcards.genome}.v6.4.zip
#             unzip -d {params.out_dir} {params.out_dir}.zip
#             mv {params.out_dir}/data/genomes/{wildcards.genome}/{wildcards.genome}.full.annotation {output}
#             rm -rf {params.out_dir} && rm -rf {params.out_dir}.zip) &> {params.check}
#         """
#
#
#
#
#     if paired_end:
#         bwa_reads=["{outdir}/fastqs/{sample}_{exp_type}_R1.fastq", "{outdir}/fastqs/{sample}_{exp_type}_R2.fastq"]
#     else:
#         bwa_reads=["{outdir}/fastqs/{sample}_{exp_type}.fastq"]
#
#
#
#     rule bam_stats:
#         input: "{outdir}/tmp/bams/{sample}.{genome}.{exp_type}.bam"
#         output: "{outdir}/info/{sample}.{genome}.{exp_type}.bam_stats.txt"
#         threads: threads
#         log: "{outdir}/logs/{sample}_{genome}_{exp_type}__bam_stats.log"
#         shell: "(samtools flagstat -@ {threads} {input} > {output}) &> {log}"
#
# if sample_type != "bigWig" and sample_type != "bedGraph":
#
#     if sample_type == "bam":
#         rule copy_bam_exp:
#             input: experiments
#             output: bam_output[0]
#             shell: """
#                 cp {input} {output}
#             """
#
#         if controls != "None":
#             rule copy_bam_ctr:
#                 input: controls
#                 output: bam_output[1]
#                 shell: """
#                     cp {input} {output}
#                 """
#
#     rule prep_bam:
#         input: "{outdir}/tmp/bams/{sample}.{genome}.{exp_type}.bam"
#         output: "{outdir}/bams_unstranded/{sample}.{genome}.{exp_type}.bam"
#         threads: threads
#         params:
#             downsample=downsample
#         log: "{outdir}/logs/{sample}_{genome}_{exp_type}__prep_bam.log"
#         shell: """
#             (
#             downsample={params.downsample}
#             if [[ "$downsample" -gt 0 ]]; then
#                 echo "Downsampling to "$downsample" reads."
#                 numread=$(samtools view -@ {threads} {input} | wc -l)
#                 echo {params.downsample}
#                 echo $numread
#                 if [[ $downsample -gt $numread ]] || [[ $downsample -eq $numread ]]; then
#                     echo "Attempted to downsample the same number or more reads than available in .bam. All reads will be used!"
#                     cp {input} {output}
#                 else
#                     downfrac=$(echo "scale=20; $downsample/$numread" | bc)
#                     echo "Downsample fraction is "$downfrac
#                     picard-tools DownsampleSam I={input} O={output} P=$downfrac VALIDATION_STRINGENCY=SILENT R=42
#                 fi
#             else
#                 echo "No downsampling. Moving bam to output."
#                 cp {input} {output}
#             fi
#             ) &> {log}
#         """
#
#     rule index_bam:
#         input: "{outdir}/bams_unstranded/{sample}.{genome}.{exp_type}.bam"
#         output: "{outdir}/bams_unstranded/{sample}.{genome}.{exp_type}.bam.bai"
#         threads: threads
#         log: "{outdir}/logs/{sample}_{genome}_{exp_type}__index_bam.log"
#         shell: """
#             (samtools index -@ {threads} {input}) &> {log}
#         """
#
#     if strand_specific and paired_end:
#         rule split_strands_pe:
#             input:
#                 bam="{outdir}/bams_unstranded/{sample}.{genome}.{exp_type}.bam"
#             output:
#                 plus="{outdir}/bams_stranded/{sample}.{genome}.{exp_type}.p.bam",
#                 minus="{outdir}/bams_stranded/{sample}.{genome}.{exp_type}.m.bam"
#             threads: threads
#             log: "{outdir}/logs/{sample}_{genome}_{exp_type}__split_strands_pe.log"
#             shell: """
#                 # Adapted from https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
#                 (mkdir {wildcards.outdir}/bams_tmp) &> /dev/null || true
#                 (samtools view -@ {threads} -b -f 128 -F 16 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.fwd1.bam
#                 samtools view -@ {threads} -b -f 80 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.fwd2.bam
#                 samtools merge -@ {threads} -f {output.plus} {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.fwd1.bam {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.fwd2.bam
#                 samtools index -@ {threads} {output.plus}
#                 samtools view -@ {threads} -b -f 144 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.rev1.bam
#                 samtools view -@ {threads} -b -f 64 -F 16 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.rev2.bam
#                 samtools merge -@ {threads} -f {output.minus} {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.rev1.bam {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.rev2.bam
#                 samtools index -@ {threads} {output.minus}
#                 rm -rf {wildcards.outdir}/bams_tmp) &> {log}
#             """
#     elif strand_specific:
#         rule split_strands_se:
#             input: "{outdir}/bams_unstranded/{sample}.{genome}.{exp_type}.bam",
#             output:
#                 plus="{outdir}/bams_stranded/{sample}.{genome}.{exp_type}.p.bam",
#                 minus="{outdir}/bams_stranded/{sample}.{genome}.{exp_type}.m.bam"
#             threads: threads
#             log: "{outdir}/logs/{sample}_{genome}_{exp_type}__split_strands_se.log"
#             shell: """
#                 # Adapted from https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
#                 (mkdir {wildcards.outdir}/bams_tmp) &> /dev/null || true
#                 (samtools view -@ {threads} -b -f 16 {input} > {output.plus}
#                 samtools index -@ {threads} {output.plus}
#                 samtools view -@ {threads} -b -F 16 {input} > {output.minus}
#                 samtools index -@ {threads} {output.minus}
#                 rm -rf {wildcards.outdir}/bams_tmp) &> {log}
#             """
#
#     rule macs2_predictd:
#         input: "{outdir}/bams_unstranded/{sample}.{genome}.{exp_type}.bam"
#         output: "{outdir}/info/{sample}.{genome}.{exp_type}_predictd.txt"
#         params: "-g " + str(effective_genome_size) + " --outdir " + outdir + "/info"
#         log: "{outdir}/logs/{sample}_{genome}_{exp_type}__macs2_predictd.log"
#         shell: """
#         (echo {params} && macs2 predictd -i {input} --outdir {wildcards.outdir}/tmp &> {output} && \
#         if grep -q "Can't find enough pairs of symmetric peaks to build model!" {output}; \
#         then echo "Failed due to predictd! Will now attempt to rerun with broaded MFOLD. See the message below:\n" && cat {output}; fi && \
#         macs2 predictd -m 1 100 -i {input} --outdir {wildcards.outdir}/tmp &> {output} && \
#         if grep -q "Can't find enough pairs of symmetric peaks to build model!" {output}; \
#         then echo "Failed due to predictd again! See the message below:\n" && cat {output} && exit 1; fi) &> {log}
#         """
#
#     # Set up arguments for peak callers
#     treat_file_macs2="{outdir}/bams_unstranded/{sample}.{genome}.experiment.bam"
#
#     if paired_end:
#         params_macs2="-f BAMPE --broad-cutoff .01 -q .01 --broad -g " + str(effective_genome_size)
#         treat_file_epic="{outdir}/tmp/pe_beds/{sample}.{genome}.experiment.bed"
#         epic_info=treat_file_epic
#         index_epic=treat_file_epic
#     else:
#         params_macs2="--broad-cutoff .01 -q .01 --broad -g " + str(effective_genome_size)
#         epic_info="{outdir}/info/{sample}.{genome}.experiment_predictd.txt"
#         treat_file_epic="{outdir}/bams_unstranded/{sample}.{genome}.experiment.bam"
#
#     if controls != "None":
#         control_file_macs2="{outdir}/bams_unstranded/{sample}.{genome}.control.bam"
#         macs2_command="(macs2 callpeak -t {input.treatment} -c {input.control}" \
#                       + " {params} --outdir {wildcards.outdir}/peaks_macs_unstranded/ " \
#                       + "-n {wildcards.sample}_{wildcards.genome}.unstranded" \
#                       + " && mv {wildcards.outdir}/peaks_macs_unstranded/{wildcards.sample}_{wildcards.genome}.unstranded_peaks.broadPeak {wildcards.outdir}/peaks_macs_unstranded/{wildcards.sample}_{wildcards.genome}.unstranded.broadPeak) &> {log}"
#         epic2_command="epic2 -t {input.treatment} -c {input.control} {params.epic} --output {output}"
#         if paired_end:
#             control_file_epic="{outdir}/tmp/pe_beds/{sample}.{genome}.control.bed"
#         else:
#             control_file_epic="{outdir}/bams_unstranded/{sample}.{genome}.control.bam"
#             index_epic=["{outdir}/bams_unstranded/{sample}.{genome}.experiment.bam.bai",
#                         "{outdir}/bams_unstranded/{sample}.{genome}.control.bam.bai"]
#     else:
#         control_file_macs2=treat_file_macs2
#         control_file_epic=treat_file_epic
#         index_epic="{outdir}/bams_unstranded/{sample}.{genome}.experiment.bam.bai"
#         macs2_command="(macs2 callpeak -t {input.treatment} {params} --outdir {wildcards.outdir}/peaks_macs_unstranded/" \
#                       + " -n {wildcards.sample}_{wildcards.genome}.unstranded" \
#                         + " && mv {wildcards.outdir}/peaks_macs_unstranded/{wildcards.sample}_{wildcards.genome}.unstranded_peaks.broadPeak {wildcards.outdir}/peaks_macs_unstranded/{wildcards.sample}_{wildcards.genome}.unstranded.broadPeak) &> {log}"
#         epic2_command="epic2 -t {input.treatment} {params.epic} --output {output}"
#
#     if paired_end:
#         epic2_command_lt="(wget {params.wget_in} -O {params.wget_out} && " + epic2_command + ") &> {log}"
#     else:
#         epic2_command_lt="(frag_med=$(cat {input.info} | grep -o 'predicted fragment length is [0-9]* bps' | cut -d ' ' -f 5) &&" \
#             + " wget {params.wget_in} -O {params.wget_out} && " + epic2_command + " --fragment-size ${{frag_med%.*}}) &> {log}"
#
#     rule macs2_callpeak_unstranded:
#         input:
#             treatment=treat_file_macs2,
#             control=control_file_macs2
#         output:
#             "{outdir}/peaks_macs_unstranded/{sample}_{genome}.unstranded.broadPeak",
#         params: params_macs2
#         log: "{outdir}/logs/{sample}_{genome}__macs2_callpeak_unstranded.log"
#         shell: macs2_command
#
#
#     rule bampe_to_bedpe:
#         # adapted from https://github.com/biocore-ntnu/epic2/issues/24
#         # awk only accepts MAPQ score > 30 and non-negative start position
#         input:
#             "{outdir}/bams_unstranded/{sample}.{genome}.{exp_type}.bam"
#         output:
#             temp("{outdir}/tmp/pe_beds/{sample}.{genome}.{exp_type}.bed")
#         params:
#             samtools_sort_extra="-n",
#             bedtools_bamtobed_extra="-bedpe"
#         threads: threads
#         log: "{outdir}/logs/{sample}_{genome}_{exp_type}__bampe_to_bedpe.log"
#         shell: """
#             (samtools sort {params.samtools_sort_extra} -@ {threads} {input} | \
#             (bedtools bamtobed -i /dev/stdin {params.bedtools_bamtobed_extra} > /dev/stdout) 2> /dev/null | \
#             awk '{{if($2 >= 1 && $8 >= 30) print}}' /dev/stdin > {output}) &> {log}
#         """
#
#
#     rule epic2_callpeak_unstranded:
#         input:
#              treatment=treat_file_epic,
#              control=control_file_epic,
#              index_epic=index_epic,
#              info=epic_info
#         output:
#             "{outdir}/peaks_epic_unstranded/{sample}_{genome}.unstranded.bed"
#         params:
#             epic="-e 100 -fdr .01 --effective-genome-fraction " + str(effective_genome_fraction) \
#                  + " --mapq 30 --chromsizes {outdir}/tmp/" + genome + ".chrom.sizes",
#             wget_in="ftp://hgdownload.soe.ucsc.edu/goldenPath/" + genome +\
#                  "/bigZips/" + genome + ".chrom.sizes",
#             wget_out="{outdir}/tmp/" + genome + ".chrom.sizes"
#         log: "{outdir}/logs/{sample}_{genome}__epic2_callpeak_unstranded.log"
#         shell: epic2_command_lt
#
#
#     # Calculate this for all reads
#     rule deeptools_coverage_unstranded:
#             input:
#                 bam="{outdir}/bams_unstranded/{sample}.{genome}.experiment.bam",
#                 index="{outdir}/bams_unstranded/{sample}.{genome}.experiment.bam.bai"
#             output: "{outdir}/coverage_unstranded/{sample}.{genome}.bw"
#             threads: threads
#             log: "{outdir}/logs/{sample}_{genome}__deeptools_coverage_unstranded.log"
#             params:
#                 extra="--ignoreForNormalization chrX chrY chrM --ignoreDuplicates --minMappingQuality" \
#                       + " 30 --binSize 10 --effectiveGenomeSize " + str(effective_genome_size)
#             shell: """
#                 (bamCoverage -b {input.bam} -p {threads} {params.extra} -o {output}) &> {log}
#             """
#
#     # Handles stranded peak calls and coverage
#     if strand_specific:
#         experiment_plus = "{outdir}/bams_stranded/{sample}.{genome}.experiment.p.bam"
#         experiment_minus = "{outdir}/bams_stranded/{sample}.{genome}.experiment.m.bam"
#         index = ["{outdir}/bams_stranded/{sample}.{genome}.experiment.p.bam.bai",
#                  "{outdir}/bams_stranded/{sample}.{genome}.experiment.m.bam.bai"]
#         if controls != "None":
#             control_plus = "{outdir}/bams_stranded/{sample}.{genome}.control.p.bam"
#             control_minus = "{outdir}/bams_stranded/{sample}.{genome}.control.m.bam"
#             index.extend(["{outdir}/bams_stranded/{sample}.{genome}.control.p.bam.bai",
#                          "{outdir}/bams_stranded/{sample}.{genome}.control.m.bam.bai"]),
#             macs2_ss_command="macs2 callpeak -t {input.experiment_plus} -c {input.control_plus} {params.callpeak} --nomodel" \
#                              + " --extsize ${{frag_med%.*}} --outdir {wildcards.outdir}/peaks_macs_stranded/ " \
#                              + "-n {wildcards.sample}_{wildcards.genome}.plus && " \
#                              + "macs2 callpeak -t {input.experiment_minus} -c {input.control_minus} {params.callpeak} --nomodel" \
#                              + " --extsize ${{frag_med%.*}} --outdir {wildcards.outdir}/peaks_macs_stranded/ " \
#                              + "-n {wildcards.sample}_{wildcards.genome}.minus"
#             epic2_ss_command="epic2 -t {input.experiment_plus} -c {input.control_plus} {params.epic} -fs ${{frag_med%.*}}" \
#                              + " --output {output.plus} && " \
#                              + "epic2 -t {input.experiment_minus} -c {input.control_minus} {params.epic} -fs ${{frag_med%.*}}" \
#                              + " --output {output.minus}"
#         else:
#             control_plus = experiment_plus
#             control_minus = experiment_minus
#             macs2_ss_command="macs2 callpeak -t {input.experiment_plus} {params.callpeak} --nomodel" \
#                              + " --extsize ${{frag_med%.*}} --outdir {wildcards.outdir}/peaks_macs_stranded/ " \
#                              + "-n {wildcards.sample}_{wildcards.genome}.plus && " \
#                              + "macs2 callpeak -t {input.experiment_minus} {params.callpeak} --nomodel" \
#                              + " --extsize ${{frag_med%.*}} --outdir {wildcards.outdir}/peaks_macs_stranded/ " \
#                              + "-n {wildcards.sample}_{wildcards.genome}.minus"
#             epic2_ss_command="epic2 -t {input.experiment_plus} {params.epic} -fs ${{frag_med%.*}}" \
#                              + " --output {output.plus} && " \
#                              + "epic2 -t {input.experiment_minus} {params.epic} -fs ${{frag_med%.*}}" \
#                              + " --output {output.minus}"
#
#         rule deeptools_coverage_stranded:
#                 input:
#                     bam="{outdir}/bams_unstranded/{sample}.{genome}.experiment.bam",
#                     index="{outdir}/bams_unstranded/{sample}.{genome}.experiment.bam.bai"
#                 output:
#                     plus="{outdir}/coverage_stranded/{sample}.{genome}.p.bw",
#                     minus="{outdir}/coverage_stranded/{sample}.{genome}.m.bw"
#                 threads: threads
#                 log: "{outdir}/logs/{sample}_{genome}__deeptools_coverage_stranded.log"
#                 params:
#                     extra="--ignoreForNormalization chrX chrY chrM --ignoreDuplicates --minMappingQuality" \
#                           + " 30 --binSize 10 --effectiveGenomeSize " + str(effective_genome_size)
#                 shell: """
#                     (bamCoverage -b {input.bam} -p {threads} {params.extra} --filterRNAstrand forward -o {output.plus}
#                     bamCoverage -b {input.bam} -p {threads} {params.extra} --filterRNAstrand reverse -o {output.minus}) &> {log}
#                 """
#
#         if paired_end:
#             rule deeptools_get_pe_fragment_sizes:
#                 # TODO: First, remove the low MAPQ score reads as they increase the overall size
#                 input:
#                     bam="{outdir}/bams_unstranded/{sample}.{genome}.experiment.bam",
#                     index="{outdir}/bams_unstranded/{sample}.{genome}.experiment.bam.bai"
#                 output:
#                     info="{outdir}/info/{sample}.{genome}.experiment.frag_lengths.txt",
#                     plot="{outdir}/info/{sample}.{genome}.experiment.frag_lengths.png"
#                 params: "--samplesLabel " + sample_name
#                 log: "{outdir}/logs/{sample}_{genome}__deeptools_get_pe_fragment_sizes.log"
#                 threads: threads
#                 shell: "(bamPEFragmentSize --bamfiles {input.bam} --histogram {output.plot} --table {output.info} -p {threads}) &> {log}"
#
#             rule epic_callpeaks_pe_stranded:
#                 input:
#                     info="{outdir}/info/{sample}.{genome}.experiment.frag_lengths.txt",
#                     experiment_plus=experiment_plus,
#                     experiment_minus=experiment_minus,
#                     control_plus=control_plus,
#                     control_minus=control_minus
#                 output:
#                     plus=temp("{outdir}/peaks_epic_stranded/{sample}_{genome}.plus.bed"),
#                     minus=temp("{outdir}/peaks_epic_stranded/{sample}_{genome}.minus.bed")
#                 params:
#                     epic="-e 100 -fdr .01 --effective-genome-fraction " + str(effective_genome_fraction) \
#                          + " --mapq 30 --chromsizes {outdir}/tmp/" + genome + ".chrom.sizes",
#                     wget_in="ftp://hgdownload.soe.ucsc.edu/goldenPath/" + genome +\
#                          "/bigZips/" + genome + ".chrom.sizes",
#                     wget_out="{outdir}/tmp/" + genome + ".chrom.sizes"
#                 log: "{outdir}/logs/{sample}_{genome}__epic_callpeaks_pe_stranded.log"
#                 shell: "(frag_med=$(head -n 2 {input.info} | tail -n 1 | awk '{{print $6}}')" \
#                        + " && wget {params.wget_in} -O {params.wget_out} && " + epic2_ss_command + ") &> {log}"
#
#             rule macs2_callpeaks_pe_stranded:
#                 input:
#                     info="{outdir}/info/{sample}.{genome}.experiment.frag_lengths.txt",
#                     experiment_plus=experiment_plus,
#                     experiment_minus=experiment_minus,
#                     control_plus=control_plus,
#                     control_minus=control_minus
#                 output:
#                     plus=temp("{outdir}/peaks_macs_stranded/{sample}_{genome}.plus_peaks.broadPeak"),
#                     minus=temp("{outdir}/peaks_macs_stranded/{sample}_{genome}.minus_peaks.broadPeak")
#                 params:
#                     callpeak="--broad-cutoff .01 -q .01 --broad -g " + str(effective_genome_size)
#                 log: "{outdir}/logs/{sample}_{genome}__macs2_callpeaks_pe_stranded.log"
#                 shell: "(frag_med=$(head -n 2 {input.info} | tail -n 1 | awk '{{print $6}}') && " \
#                         + macs2_ss_command + ") &> {log}"
#
#         else:
#             rule epic_callpeaks_se_stranded:
#                 input:
#                     info="{outdir}/info/{sample}.{genome}.experiment_predictd.txt",
#                     experiment_plus=experiment_plus,
#                     experiment_minus=experiment_minus,
#                     control_plus=control_plus,
#                     control_minus=control_minus
#                 output:
#                     plus=temp("{outdir}/peaks_epic_stranded/{sample}_{genome}.plus.bed"),
#                     minus=temp("{outdir}/peaks_epic_stranded/{sample}_{genome}.minus.bed")
#                 params:
#                     epic="-e 100 -fdr .01 --effective-genome-fraction " + str(effective_genome_fraction) \
#                          + " --mapq 30 --chromsizes {outdir}/tmp/" + genome + ".chrom.sizes",
#                     wget_in="ftp://hgdownload.soe.ucsc.edu/goldenPath/" + genome +\
#                          "/bigZips/" + genome + ".chrom.sizes",
#                     wget_out="{outdir}/tmp/" + genome + ".chrom.sizes"
#                 log: "{outdir}/logs/{sample}_{genome}__epic_callpeaks_se_stranded.log"
#                 shell: "(frag_med=$(cat {input.info} | grep -o 'predicted fragment length is [0-9]* bps' | cut -d ' ' -f 5)" \
#                        + " && wget {params.wget_in} -O {params.wget_out} && " + epic2_ss_command + ") &> {log}"
#
#             rule macs2_callpeaks_se_stranded:
#                 # Based on https://github.com/PEHGP/drippipline
#                 input:
#                     info="{outdir}/info/{sample}.{genome}.experiment_predictd.txt",
#                     experiment_plus=experiment_plus,
#                     experiment_minus=experiment_minus,
#                     control_plus=control_plus,
#                     control_minus=control_minus
#                 output:
#                     plus=temp("{outdir}/peaks_macs_stranded/{sample}_{genome}.plus_peaks.broadPeak"),
#                     minus=temp("{outdir}/peaks_macs_stranded/{sample}_{genome}.minus_peaks.broadPeak")
#                 params:
#                     callpeak="--broad-cutoff .01 -q .01 --broad -g " + str(effective_genome_size)
#                 log: "{outdir}/logs/{sample}_{genome}__macs2_callpeaks_se_stranded.log"
#                 shell: "(frag_med=$(cat {input.info} | grep -o 'predicted fragment length is [0-9]* bps' | cut -d ' ' -f 5)" \
#                        + " && " + macs2_ss_command + ") &> {log}"
#
#         rule merge_stranded_peaks:
#             input:
#                 macs2_plus="{outdir}/peaks_macs_stranded/{sample}_{genome}.plus_peaks.broadPeak",
#                 macs2_minus="{outdir}/peaks_macs_stranded/{sample}_{genome}.minus_peaks.broadPeak",
#                 epic2_plus="{outdir}/peaks_epic_stranded/{sample}_{genome}.plus.bed",
#                 epic2_minus="{outdir}/peaks_epic_stranded/{sample}_{genome}.minus.bed"
#             output:
#                 macs2_sorted="{outdir}/peaks_macs_stranded/{sample}_{genome}.stranded.broadPeak",
#                 epic2_sorted="{outdir}/peaks_epic_stranded/{sample}_{genome}.stranded.bed"
#             params:
#                 epic2_unsorted="{outdir}/peaks_epic_stranded/{sample}_{genome}.stranded.unsorted.bed",
#                 macs2_unsorted="{outdir}/peaks_macs_stranded/{sample}_{genome}.stranded.unsorted.broadPeak"
#             shell: """
#             awk ' {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""+""\t"$7"\t"$8"\t"$9}} ' {input.macs2_plus} > {params.macs2_unsorted}
#             awk ' {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""-""\t"$7"\t"$8"\t"$9}} ' {input.macs2_minus} >> {params.macs2_unsorted}
#             bedtools sort -i {params.macs2_unsorted} > {output.macs2_sorted}
#             awk ' NR==1 {{print; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""+""\t"$7"\t"$8"\t"$9"\t"$10}} ' {input.epic2_plus} > {params.epic2_unsorted}
#             awk ' NR==1 {{print; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""-""\t"$7"\t"$8"\t"$9"\t"$10}} ' {input.epic2_minus} >> {params.epic2_unsorted}
#             bedtools sort -i {params.epic2_unsorted} > {output.epic2_sorted}
#             """
#
# ## TODO: Correlation module goes here
# rule get_correlation_bin_sthreads:
#     input: "{outdir}/coverage_unstranded/{sample}.{genome}.bw",
#     params:
#         helpers_dir = helpers_dir
#     output:
#         npz=temp("{outdir}/QC/{sample}.{genome}.gold_standard_bin_sthreads.npz"),
#         tab="{outdir}/QC/{sample}.{genome}.gold_standard_bin_sthreads.tab"
#     threads: threads
#     log: "{outdir}/logs/{sample}_{genome}__get_correlation_bin_sthreads.log"
#     shell:
#         """
#         (multiBigwigSummary BED-file --BED {params.helpers_dir}/data/correlation_genes_100kb.{wildcards.genome}.1kbwindow.bed \
#         -o {output.npz} -b {input} --outRawCounts {output.tab} -p {threads}) &> {log}
#         """
#
#
# rule correlation_analysis:
#     input: "{outdir}/QC/{sample}.{genome}.gold_standard_bin_sthreads.tab"
#     output:
#         image="{outdir}/QC/{sample}.{genome}.correlation.png",
#         data="{outdir}/QC/{sample}.{genome}.correlation.rda"
#     params:
#         helpers_dir = helpers_dir,
#         mode = mode
#     log: "{outdir}/logs/{sample}_{genome}__correlation_analysis.log"
#     shell:
#          """
#          (Rscript {params.helpers_dir}/correlation_test.R {input} {wildcards.sample} {params.mode} {params.helpers_dir}) &> {log}
#          """
#
#
# rule assign_genome_annotations:
#     # Calculates the percentage of called peaks overlapping with repeat regions
#     input:
#         gene_anno=genome_home_dir + "/{genome}/homer_anno.txt",
#         peaks="{outdir}/peaks_final_{strand}/{sample}_{genome}.{strand}.bed"
#     output:
#         stats_out="{outdir}/QC/{sample}_{genome}_{strand}.feature_overlaps.txt"
#     log: "{outdir}/logs/{sample}_{genome}_{strand}__assign_genome_annotations.log"
#     shell: """
#         (assignGenomeAnnotation {input.peaks} {input.gene_anno} > {output.stats_out}) &> {log}
#     """
#
#
# rule prepare_report:
#     input: final_report_input
#     output:
#         html="{outdir}/{sample}_{genome}.QC_report.html",
#         rda="{outdir}/{sample}_{genome}.QC_report.rda"
#     params:
#         helpers_dir=helpers_dir,
#         sample_name=sample_name,
#         configs=outdirorig + "/rseqVars.json",
#         final_report_dict_file=final_report_dict_file[0]
#     log: "{outdir}/logs/{sample}_{genome}__prepare_report.log"
#     shell:
#      """
#      (echo {output}
#      Rscript {params.helpers_dir}/prepare_report.R {params.final_report_dict_file} {params.sample_name} {params.configs}
#      rm {params.final_report_dict_file}) &> {log}
#      """
#
#
# # rule run_QmRLFS_finder:
# #     # Packaged with RSeq from: https://github.com/piroonj/QmRLFS-finder (Oct 12 2020)
# #     input: genome_home_dir + "/{genome}/{genome}.fa"
# #     output:
# #         table=genome_home_dir + "/{genome}/rloop_predictions/RLFS.{genome}.out.table.txt",
# #         bed=genome_home_dir + "/{genome}/rloop_predictions/RLFS.{genome}.out.table.bed"
# #     params:
# #         helpers_dir=helpers_dir,
# #         bed=genome_home_dir + "/{genome}/rloop_predictions/RLFS.{genome}.out.table.raw.bed",
# #         genome_home_dir=genome_home_dir
# #     shell: """
# #     python {params.helpers_dir}/external/QmRLFS-finder.py -i {input} \
# #     -o {params.genome_home_dir}/{wildcards.genome}/rloop_predictions/RLFS.{wildcards.genome}
# #     awk '{{OFS="\t";print($1,$4,$14,$3,0,$21)}}' {output.table} > {params.bed}
# #     bedtools sort -i {params.bed} | mergeBed -i stdin -s | awk '{{OFS="\t";print($1,$2,$3,".",".",$4)}}' > {output.bed}
# #     """
#
# rule download_rlfs_annotations:
#     output:
#         bed=genome_home_dir + "/{genome}/rloop_predictions/{genome}.rlfs.bed"
#     shell: """
#     wget -O {output} https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/{wildcards.genome}.rlfs.bed
#     """
#
#
# rule rlfs_enrichment:
#     input:
#          peaks="{outdir}/peaks_final_unstranded/{sample}_{genome}.unstranded.bed",
#          rlfs=genome_home_dir + "/{genome}/rloop_predictions/{genome}.rlfs.bed"
#     output: "{outdir}/QC/{sample}_{genome}.rlfs_enrichment.rda"
#     params:
#         helpers_dir=helpers_dir
#     threads: threads
#     log: "{outdir}/logs/{sample}_{genome}__rlfs_enrichment.log"
#     shell: """
#      (Rscript {params.helpers_dir}/rlfs_perm_test.R {threads} {wildcards.genome} {input.peaks} {input.rlfs} {output}) &> {log}
#      """
#
#
# rule compile_peaks:
#     input:
#         macs2="{outdir}/peaks_macs_{strand}/{sample}_{genome}.{strand}.broadPeak",
#         epic2="{outdir}/peaks_epic_{strand}/{sample}_{genome}.{strand}.bed"
#     output:
#         output_bed="{outdir}/peaks_final_{strand}/{sample}_{genome}.{strand}.bed",
#         output_rda="{outdir}/peaks_final_{strand}/{sample}_{genome}.{strand}.rda"
#     params:
#         configs=outdirorig + "/rseqVars.json",
#         output_bed_tmp="{outdir}/peaks_final_{strand}/{sample}_{genome}.{strand}.tmp.bed",
#         sample_name=sample_name,
#         helpers_dir=helpers_dir
#     log: "{outdir}/logs/{sample}_{genome}_{strand}__compile_peaks_{strand}.log"
#     shell: """
#     (
#     Rscript {params.helpers_dir}/compile_peaks.R {params.configs} {params.sample_name} {input.macs2} {input.epic2} {output.output_rda} {params.output_bed_tmp}
#     awk ' $2 >= 0 ' {params.output_bed_tmp} > {output.output_bed} && rm {params.output_bed_tmp}
#     ) &> {log}
#     """
#
#
#
#
#
#
