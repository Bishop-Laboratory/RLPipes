########################################################################################################################
##############################################   Parse inputs    #######################################################
########################################################################################################################

import json
import os

# Config vars
helpers_dir=config['helpers_dir'][0]
mode=config['mode'][0]
paired_end=config['paired_end'][0]
strand_specific=config['strand_specific'][0]
homer_anno_available=config['homer_anno_available'][0]
genome=config['genome'][0]
genome_home_dir=config['genome_home_dir'][0]
effective_genome_size=config['effective_genome_size'][0]
full_genome_length=config['full_genome_length'][0]
effective_genome_fraction = effective_genome_size/full_genome_length
cores=int(config['cores'][0])
sample_type=config['sample_type'][0]
no_dedupe=config['no_dedupe'][0]
no_fastp=config['no_fastp'][0]
sample_name=config['sample_name'][0]
# TODO: Should there be an option to specify whether or not to separate output folders based on sample_name?
outdirorig=config['out_dir'][0]
outdir=config['out_dir'][0] + "/" + sample_name
outdir=outdir.replace("//", "/")
experiments=config['experiments'][0].split(",")[0]
controls=config['controls'][0]
if controls != "None":
    controls=controls.split(",")[0]

# Set expected output and inputs for merging of technical replicates (assumes public samples)
output_fastq_experiment_1 = expand("{outdir}/fastqs/{sample_name}_experiment_R1.fastq",
                 sample_name=sample_name, outdir=outdir)
output_fastq_experiment_2 = expand("{outdir}/fastqs/{sample_name}_experiment_R2.fastq",
                 sample_name=sample_name, outdir=outdir)
if paired_end:
    if sample_type == "public":
        merge_input_experiment_1 = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}_1.fastq",
                                   sample_name=sample_name,
                                   srr_acc=experiments, outdir=outdir)
        merge_input_experiment_2 = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}_2.fastq",
                                   sample_name=sample_name,
                                   srr_acc=experiments, outdir=outdir)
        merge_input_control_1 = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}_1.fastq",
                               sample_name=sample_name,
                               srr_acc=controls, outdir=outdir)
        merge_input_control_2 = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}_2.fastq",
                                   sample_name=sample_name,
                                   srr_acc=controls, outdir=outdir)
    elif sample_type == "fastq":
        # Fastq files are supplied, and they are paired end
        merge_input_experiment_1 = experiments.split("+")[0]
        merge_input_experiment_2 = experiments.split("+")[1]
        if controls == "None":
            merge_input_control_1 = None
            merge_input_control_2 = None
        else:
            merge_input_control_1 = controls.split("+")[0]
            merge_input_control_2 = controls.split("+")[1]
else:
    if sample_type == "public":
        merge_input_experiment = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}_1.fastq",
                                   sample_name=sample_name,
                                   srr_acc=experiments, outdir=outdir)
        merge_input_control = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}_1.fastq",
                               sample_name=sample_name,
                               srr_acc=controls, outdir=outdir)
    else:
        # Fastq files are supplied, and they are single end
        merge_input_experiment = experiments
        if controls == "None":
            merge_input_control = None
        else:
            merge_input_control = controls


output_fastq_control_1 = expand("{outdir}/fastqs/{sample_name}_control_R1.fastq",
                 sample_name=sample_name, outdir=outdir)
output_fastq_control_2 = expand("{outdir}/fastqs/{sample_name}_control_R2.fastq",
                 sample_name=sample_name, outdir=outdir)

output_fastq_experiment = expand("{outdir}/fastqs/{sample_name}.fastq",
                 sample_name=sample_name, outdir=outdir)
output_fastq_control = expand("{outdir}/fastqs/{sample_name}.fastq",
                 sample_name=sample_name, outdir=outdir)

# Set expected outputs of runs
if strand_specific:
    coverage_output = [expand("{outdir}/coverage_stranded/{sample}.{genome}.p.bw", genome=genome,
                        sample=sample_name, outdir=outdir),
                        expand("{outdir}/coverage_unstranded/{sample}.{genome}.bw", genome=genome,
                            sample=sample_name, outdir=outdir)]
    peak_output = [expand("{outdir}/peaks_macs_unstranded/{sample}_{genome}_peaks.broadPeak", genome=genome,
                         sample=sample_name, outdir=outdir),
                   expand("{outdir}/peaks_epic_unstranded/{sample}_{genome}.bed", genome=genome,
                         sample=sample_name, outdir=outdir),
                   expand("{outdir}/peaks_epic_stranded/{sample}_{genome}_plus.bed", genome=genome,
                         sample=sample_name, outdir=outdir),
                   expand("{outdir}/peaks_epic_stranded/{sample}_{genome}_minus.bed", genome=genome,
                         sample=sample_name, outdir=outdir),
                   expand("{outdir}/peaks_macs_stranded/{sample}_{genome}_plus_peaks.xls", genome=genome,
                         sample=sample_name, outdir=outdir)]
else:
    coverage_output = expand("{outdir}/coverage_unstranded/{sample}.{genome}.bw", genome=genome,
                            sample=sample_name, outdir=outdir)
    peak_output = [expand("{outdir}/peaks_macs_unstranded/{sample}_{genome}_peaks.broadPeak", genome=genome,
                         sample=sample_name, outdir=outdir),
                   expand("{outdir}/peaks_epic_unstranded/{sample}_{genome}.bed", genome=genome,
                         sample=sample_name, outdir=outdir)
                   ]

if genome == "hg38":
    correlation_output = expand("{outdir}/QC/{sample}.{genome}.correlation.{file_ends}", genome=genome,
                                sample=sample_name, outdir=outdir, file_ends=["png", "rda"])
    rlfs_output = expand("{outdir}/QC/{sample}_{genome}.rlfs_enrichment.rda", genome=genome,
                                sample=sample_name, outdir=outdir)
else:
    correlation_output = peak_output
    rlfs_output = peak_output

if homer_anno_available:
    anno_output = expand("{outdir}/QC/{sample}_{genome}.feature_overlaps.txt", genome=genome,
                                sample=sample_name, outdir=outdir)
    # rlfs_output = expand("{outdir}/QC/{sample}_{genome}.rlfs_enrichment.rda", genome=genome,
    #                             sample=sample_name, outdir=outdir)
else:
    anno_output = peak_output
    # rlfs_output = peak_output

if controls != "None":
    bam_output = expand("{outdir}/bams/{sample_name}.{genome}.{exp_type}.bam",
                            sample_name=sample_name, outdir=outdir, genome=genome, exp_type=["experiment","control"])
    bam_index_output = expand("{outdir}/bams/{sample_name}.{genome}.{exp_type}.bam.bai",
                            sample_name=sample_name, outdir=outdir, genome=genome, exp_type=["experiment","control"])
    bam_stats_output = expand("{outdir}/info/{sample_name}.{genome}.{exp_type}.bam_stats.txt",
                            sample_name=sample_name, outdir=outdir, genome=genome, exp_type=["experiment","control"])
else:
    bam_output = expand("{outdir}/bams/{sample_name}.{genome}.{exp_type}.bam",
                            sample_name=sample_name, outdir=outdir, genome=genome, exp_type=["experiment",])
    bam_index_output = expand("{outdir}/bams/{sample_name}.{genome}.{exp_type}.bam.bai",
                            sample_name=sample_name, outdir=outdir, genome=genome, exp_type=["experiment",])
    bam_stats_output = expand("{outdir}/info/{sample_name}.{genome}.{exp_type}.bam_stats.txt",
                            sample_name=sample_name, outdir=outdir, genome=genome, exp_type=["experiment",])


# TODO: Put something real here
rlcons_output=peak_output

# Input for final report
final_report_input = [
    anno_output,
    correlation_output,
    expand("{outdir}/QC/fastq/{sample}.experiment.json", outdir=outdir, sample=sample_name),
    rlfs_output,
    rlcons_output,
    bam_stats_output
]
final_report_dict = {
    "anno_output": anno_output,
    "correlation_output": correlation_output,
    "fastpdata": expand("{outdir}/QC/fastq/{sample}.experiment.json", outdir=outdir, sample=sample_name),
    "rlfs_output": rlfs_output,
    "rlcons_output": rlcons_output,
    "bam_stats_output": bam_stats_output
}

os.makedirs(outdir, exist_ok=True)
final_report_dict_file = expand("{outdir}/{sample}.final_report.tmp.json", outdir=outdir, sample=sample_name)
with open(final_report_dict_file[0], "w") as write_file:
    json.dump(final_report_dict, write_file)

final_report_output = expand("{outdir}/{sample}_{genome}.QC_report.{ext}", genome=genome,
                                sample=sample_name, outdir=outdir, ext=['html', 'rda'])



# TODO: Set the TMPDIR

########################################################################################################################
##############################################   Main pipeline    ######################################################
########################################################################################################################
rule output:
    input:
        peak_output,
        coverage_output,
        bam_output,
        bam_index_output,
        bam_stats_output,
        correlation_output,
        anno_output,
        final_report_output


if sample_type != "bam" and sample_type != "bigWig" and sample_type != "bedGraph":
    if sample_type == "public":
        rule download_sra:
            threads:cores,
            output: temp("{outdir}/tmp/sras/{sample}/{srr_acc}.sra")
            log: "{outdir}/logs/{sample}_{srr_acc}__download_sra.log"
            shell: "(prefetch {wildcards.srr_acc} --output-file {output}) &> {log}"

    if paired_end:
        if sample_type == "public":
            # TODO: fastq merging only works with public samples at the moment. This should also work with local samples.
            rule sra_to_fastq_pe:
                input:
                    "{outdir}/tmp/sras/{sample}/{srr_acc}.sra"
                output:
                    temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}_1.fastq"),
                    temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}_2.fastq")
                threads: cores
                log: "{outdir}/logs/{sample}_{srr_acc}__sra_to_fastq_pe.log"
                shell:
                    "(parallel-fastq-dump -t {threads} --split-files -O" \
                    + " {wildcards.outdir}/tmp/fastqs_raw/{wildcards.sample}/ -s {input}) &> {log}"

        rule merge_fastq_pe_experiment:
            input:
                R1=merge_input_experiment_1,
                R2=merge_input_experiment_2
            output:
                R1=temp("{outdir}/tmp/fastqs_dup/{sample}_experiment_R1.fastq"),
                R2=temp("{outdir}/tmp/fastqs_dup/{sample}_experiment_R2.fastq")
            log: "{outdir}/logs/{sample}_merge_fastq_pe_experiment.log"
            shell: """
                ( cat {input.R1} > {output.R1}
                cat {input.R2} > {output.R2} ) &> {log}
            """

        if controls != "None":
            rule merge_fastq_pe_control:
                input:
                    R1=merge_input_control_1,
                    R2=merge_input_control_2
                output:
                    R1=temp("{outdir}/tmp/fastqs_dup/{sample}_control_R1.fastq"),
                    R2=temp("{outdir}/tmp/fastqs_dup/{sample}_control_R2.fastq")
                log: "{outdir}/logs/{sample}_merge_fastq_pe_control.log"
                shell: """
                    ( cat {input.R1} > {output.R1}
                    cat {input.R2} > {output.R2} ) &> {log}
                """

        if not no_fastp:
            rule fastp_pe:
                input:
                    sample=["{outdir}/tmp/fastqs_dup/{sample}_{exp_type}_R1.fastq",
                            "{outdir}/tmp/fastqs_dup/{sample}_{exp_type}_R2.fastq"]
                output:
                    trimmed=[temp("{outdir}/fastqs/{sample}_{exp_type}_R1.fastq"),
                             temp("{outdir}/fastqs/{sample}_{exp_type}_R2.fastq")],
                    html="{outdir}/QC/fastq/{sample}.{exp_type}.html",
                    json="{outdir}/QC/fastq/{sample}.{exp_type}.json"
                log: "{outdir}/logs/{sample}_{exp_type}__fastp_pe.log"
                params:
                    extra=""
                threads: cores
                wrapper:
                    "0.63.0/bio/fastp"

    else:
        if sample_type == "public":
            rule sra_to_fastq_se:
                input: "{outdir}/tmp/sras/{sample}/{srr_acc}.sra"
                output: temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}_1.fastq")
                threads: cores
                log: "{outdir}/logs/{sample}_{srr_acc}__sra_to_fastq_se.log"
                shell:
                    "(parallel-fastq-dump -t {threads} --split-files -O" \
                    + " {wildcards.outdir}/tmp/fastqs_raw/{wildcards.sample}/ -s {input}) &> {log}"

        rule merge_fastq_se_experiment:
            input: merge_input_experiment
            output: temp("{outdir}/tmp/fastqs_dup/{sample}_experiment.fastq")
            log: "{outdir}/logs/{sample}_merge_fastq_se_experiment.log"
            shell: "(cat {input} > {output}) &> {log}"

        if controls != "None":
            rule merge_fastq_se_control:
                input: merge_input_control
                output: temp("{outdir}/tmp/fastqs_dup/{sample}_control.fastq")
                log: "{outdir}/logs/{sample}_merge_fastq_se_control.log"
                shell: "(cat {input} > {output}) &> {log}"

        if not no_fastp:
            rule fastp_se:
                input:
                     sample=["{outdir}/tmp/fastqs_dup/{sample}_{exp_type}.fastq"]
                output:
                    trimmed=temp("{outdir}/fastqs/{sample}_{exp_type}.fastq"),
                    html="{outdir}/QC/fastq/{sample}.{exp_type}.html",
                    json="{outdir}/QC/fastq/{sample}.{exp_type}.json"
                log: "{outdir}/logs/{sample}_{exp_type}__fastp_se.log"
                params:
                    extra=""
                threads: cores
                wrapper:
                    "0.63.0/bio/fastp"

    rule download_fasta:
        output:
            genome_home_dir + "/{genome}/{genome}.fa"
        params:
              prefix=genome_home_dir + "/{genome}/bwa_index/{genome}",
              check=outdir + "/logs/" + sample_name + "_{genome}__download_fasta.log"
        shell: """
            (mkdir -p {params.prefix}
            wget -O {output}.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz
            gunzip {output}.gz && echo "Indexing BAM at {params.prefix} ... This may take some time ...") &> {params.check}
        """

    rule download_homer_anno:
        output:
            gene_anno=genome_home_dir + "/{genome}/homer_anno.txt"
        params:
            out_dir=genome_home_dir + "/{genome}/homer_anno",
            check=outdir + "/logs/" + sample_name + "_{genome}__download_homer_anno.log"
        shell: """
            (wget -O {params.out_dir}.zip http://homer.ucsd.edu/homer/data/genomes/{wildcards.genome}.v6.4.zip
            unzip -d {params.out_dir} {params.out_dir}.zip
            mv {params.out_dir}/data/genomes/{wildcards.genome}/{wildcards.genome}.full.annotation {output}
            rm -rf {params.out_dir} && rm -rf {params.out_dir}.zip) &> {params.check}
        """

    rule bwa_index:
        input:
            genome_home_dir + "/{genome}/{genome}.fa"
        output:
             genome_home_dir + "/{genome}/bwa_index/{genome}.amb",
             genome_home_dir + "/{genome}/bwa_index/{genome}.ann",
             genome_home_dir + "/{genome}/bwa_index/{genome}.bwt",
             genome_home_dir + "/{genome}/bwa_index/{genome}.pac",
             genome_home_dir + "/{genome}/bwa_index/{genome}.sa",
        params:
            prefix=genome_home_dir + "/{genome}/bwa_index/{genome}",
        log:
            genome_home_dir + "/{genome}/bwa_index/{genome}_bwa_index.log"
        wrapper:
            "0.63.0/bio/bwa/index"


    if paired_end:
        bwa_reads=["{outdir}/fastqs/{sample}_{exp_type}_R1.fastq", "{outdir}/fastqs/{sample}_{exp_type}_R2.fastq"]
    else:
        bwa_reads=["{outdir}/fastqs/{sample}_{exp_type}.fastq"]

    rule bwa_mem:
        input:
            bwa_index_done=genome_home_dir + "/{genome}/bwa_index/{genome}.amb",
            reads=bwa_reads
        output:
            bam="{outdir}/bams/{sample}.{genome}.{exp_type}.bam"
        log: "{outdir}/logs/{sample}_{genome}_{exp_type}__bwa_mem.log"
        params:
            index=genome_home_dir + "/{genome}/bwa_index/{genome}",
            bwa_extra=r"-R '@RG\tID:{sample}_{exp_type}\tSM:{sample}_{exp_type}'",
            samtools_sort_extra="-O BAM"
        threads: cores
        shell: """
            (bwa mem -t {threads} {params.bwa_extra} {params.index} {input.reads} | \
            samtools view -q 10 -b -@ {threads} - | \
            samtools collate -O -@ {threads} - - | \
            samtools fixmate -m - - | \
            samtools sort {params.samtools_sort_extra} -@ {threads} - | \
            samtools markdup -s -@ {threads} - {output.bam}) &> {log}
         """

    rule bam_stats:
        input: "{outdir}/bams/{sample}.{genome}.{exp_type}.bam"
        output: "{outdir}/info/{sample}.{genome}.{exp_type}.bam_stats.txt"
        threads: cores
        log: "{outdir}/logs/{sample}_{genome}_{exp_type}__bam_stats.log"
        shell: "(samtools flagstat -@ {threads} {input} > {output}) &> {log}"

if sample_type != "bigWig" and sample_type != "bedGraph":

    if sample_type == "bam":
        rule copy_bam_exp:
            input: experiments
            output: bam_output[0]
            shell: """
                cp {input} {output}
            """

        if controls != "None":
            rule copy_bam_ctr:
                input: controls
                output: bam_output[1]
                shell: """
                    cp {input} {output}
                """

    rule index_bam:
        input: "{outdir}/bams/{sample}.{genome}.{exp_type}.bam"
        output: "{outdir}/bams/{sample}.{genome}.{exp_type}.bam.bai"
        threads: cores
        log: "{outdir}/logs/{sample}_{genome}_{exp_type}__index_bam.log"
        shell: """
            (samtools index -@ {threads} {input}) &> {log}
        """

    if strand_specific and paired_end:
        rule split_strands_pe:
            input:
                bam="{outdir}/bams/{sample}.{genome}.{exp_type}.bam"
            output:
                plus="{outdir}/bams_stranded/{sample}.{genome}.{exp_type}.p.bam",
                minus="{outdir}/bams_stranded/{sample}.{genome}.{exp_type}.m.bam"
            threads: cores
            log: "{outdir}/logs/{sample}_{genome}_{exp_type}__split_strands_pe.log"
            shell: """
                # Adapted from https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
                (mkdir {wildcards.outdir}/bams_tmp) &> /dev/null || true
                (samtools view -@ {threads} -b -f 128 -F 16 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.fwd1.bam
                samtools view -@ {threads} -b -f 80 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.fwd2.bam
                samtools merge -@ {threads} -f {output.plus} {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.fwd1.bam {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.fwd2.bam
                samtools index -@ {threads} {output.plus}
                samtools view -@ {threads} -b -f 144 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.rev1.bam
                samtools view -@ {threads} -b -f 64 -F 16 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.rev2.bam
                samtools merge -@ {threads} -f {output.minus} {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.rev1.bam {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.{wildcards.exp_type}.rev2.bam
                samtools index -@ {threads} {output.minus}) &> {log}
            """
    elif strand_specific:
        rule split_strands_se:
            input: "{outdir}/bams/{sample}.{genome}.{exp_type}.bam",
            output:
                plus="{outdir}/bams_stranded/{sample}.{genome}.{exp_type}.p.bam",
                minus="{outdir}/bams_stranded/{sample}.{genome}.{exp_type}.m.bam"
            threads: cores
            log: "{outdir}/logs/{sample}_{genome}_{exp_type}__split_strands_se.log"
            shell: """
                # Adapted from https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
                (mkdir {wildcards.outdir}/bams_tmp) &> /dev/null || true
                (samtools view -@ {threads} -b -f 16 {input} > {output.plus}
                samtools index -@ {threads} {output.plus}
                samtools view -@ {threads} -b -F 16 {input} > {output.minus}
                samtools index -@ {threads} {output.minus}) &> {log}
            """

    rule macs2_predictd:
        input: "{outdir}/bams/{sample}.{genome}.{exp_type}.bam"
        output: "{outdir}/info/{sample}.{genome}.{exp_type}_predictd.txt"
        params: "-g " + str(effective_genome_size) + " --outdir " + outdir + "/info"
        log: "{outdir}/logs/{sample}_{genome}_{exp_type}__macs2_predictd.log"
        shell: """
        (echo {params} && macs2 predictd -i {input} --outdir {wildcards.outdir}/tmp &> {output} && \
        if grep -q "Can't find enough pairs of symmetric peaks to build model!" {output}; \
        then echo "Failed due to predictd! Will now attempt to rerun with broaded MFOLD. See the message below:\n" && cat {output}; fi && \
        macs2 predictd -m 1 100 -i {input} --outdir {wildcards.outdir}/tmp &> {output} && \
        if grep -q "Can't find enough pairs of symmetric peaks to build model!" {output}; \
        then echo "Failed due to predictd again! See the message below:\n" && cat {output} && exit 1; fi) &> {log}
        """

    # Set up arguments for peak callers
    treat_file_macs2="{outdir}/bams/{sample}.{genome}.experiment.bam"

    if paired_end:
        params_macs2="-f BAMPE --broad-cutoff .01 -q .01 --broad -g " + str(effective_genome_size)
        treat_file_epic="{outdir}/tmp/pe_beds/{sample}.{genome}.experiment.bed"
        epic_info=treat_file_epic
        index_epic=treat_file_epic
    else:
        params_macs2="--broad-cutoff .01 -q .01 --broad -g " + str(effective_genome_size)
        epic_info="{outdir}/info/{sample}.{genome}.experiment_predictd.txt"
        treat_file_epic="{outdir}/bams/{sample}.{genome}.experiment.bam"

    if controls != "None":
        control_file_macs2="{outdir}/bams/{sample}.{genome}.control.bam"
        macs2_command="(macs2 callpeak -t {input.treatment} -c {input.control}" \
                      + " {params} --outdir {wildcards.outdir}/peaks_macs_unstranded/ " \
                      + "-n {wildcards.sample}_{wildcards.genome}) &> {log}"
        epic2_command="epic2 -t {input.treatment} -c {input.control} {params.epic} --output {output}"
        if paired_end:
            control_file_epic="{outdir}/tmp/pe_beds/{sample}.{genome}.control.bed"
        else:
            control_file_epic="{outdir}/bams/{sample}.{genome}.control.bam"
            index_epic=["{outdir}/bams/{sample}.{genome}.experiment.bam.bai",
                        "{outdir}/bams/{sample}.{genome}.control.bam.bai"]
    else:
        control_file_macs2=treat_file_macs2
        control_file_epic=treat_file_epic
        index_epic="{outdir}/bams/{sample}.{genome}.experiment.bam.bai"
        macs2_command="(macs2 callpeak -t {input.treatment} {params} --outdir {wildcards.outdir}/peaks_macs_unstranded/" \
                      + " -n {wildcards.sample}_{wildcards.genome}) &> {log}"
        epic2_command="epic2 -t {input.treatment} {params.epic} --output {output}"

    if paired_end:
        epic2_command_lt="(wget {params.wget_in} -O {params.wget_out} && " + epic2_command + ") &> {log}"
    else:
        epic2_command_lt="(frag_med=$(cat {input.info} | grep -o 'predicted fragment length is [0-9]* bps' | cut -d ' ' -f 5) &&" \
            + " wget {params.wget_in} -O {params.wget_out} && " + epic2_command + " --fragment-size ${{frag_med%.*}}) &> {log}"

    rule macs2_callpeak_unstranded:
        input:
            treatment=treat_file_macs2,
            control=control_file_macs2
        output:
            "{outdir}/peaks_macs_unstranded/{sample}_{genome}_peaks.broadPeak",
        params: params_macs2
        log: "{outdir}/logs/{sample}_{genome}__macs2_callpeak_unstranded.log"
        shell: macs2_command


    rule bampe_to_bedpe:
        # adapted from https://github.com/biocore-ntnu/epic2/issues/24
        # awk only accepts MAPQ score > 30 and non-negative start position
        input:
            "{outdir}/bams/{sample}.{genome}.{exp_type}.bam"
        output:
            temp("{outdir}/tmp/pe_beds/{sample}.{genome}.{exp_type}.bed")
        params:
            samtools_sort_extra="-n",
            bedtools_bamtobed_extra="-bedpe"
        threads: cores
        log: "{outdir}/logs/{sample}_{genome}_{exp_type}__bampe_to_bedpe.log"
        shell: """
            (samtools sort {params.samtools_sort_extra} -@ {threads} {input} | \
            (bedtools bamtobed -i /dev/stdin {params.bedtools_bamtobed_extra} > /dev/stdout) 2> /dev/null | \
            awk '{{if($2 >= 1 && $8 >= 30) print}}' /dev/stdin > {output}) &> {log}
        """


    rule epic2_callpeak_unstranded:
        input:
             treatment=treat_file_epic,
             control=control_file_epic,
             index_epic=index_epic,
             info=epic_info
        output:
            "{outdir}/peaks_epic_unstranded/{sample}_{genome}.bed"
        params:
            epic="-e 100 -fdr .01 --effective-genome-fraction " + str(effective_genome_fraction) \
                 + " --mapq 30 --chromsizes {outdir}/tmp/" + genome + ".chrom.sizes",
            wget_in="ftp://hgdownload.soe.ucsc.edu/goldenPath/" + genome +\
                 "/bigZips/" + genome + ".chrom.sizes",
            wget_out="{outdir}/tmp/" + genome + ".chrom.sizes"
        log: "{outdir}/logs/{sample}_{genome}__epic2_callpeak_unstranded.log"
        shell: epic2_command_lt


    # Calculate this for all reads
    rule deeptools_coverage_unstranded:
            input:
                bam="{outdir}/bams/{sample}.{genome}.experiment.bam",
                index="{outdir}/bams/{sample}.{genome}.experiment.bam.bai"
            output: "{outdir}/coverage_unstranded/{sample}.{genome}.bw"
            threads: cores
            log: "{outdir}/logs/{sample}_{genome}__deeptools_coverage_unstranded.log"
            params:
                extra="--ignoreForNormalization chrX chrY chrM --ignoreDuplicates --minMappingQuality" \
                      + " 30 --binSize 10 --effectiveGenomeSize " + str(effective_genome_size)
            shell: """
                (bamCoverage -b {input.bam} -p {threads} {params.extra} -o {output}) &> {log}
            """

    # Handles stranded peak calls and coverage
    if strand_specific:
        experiment_plus = "{outdir}/bams_stranded/{sample}.{genome}.experiment.p.bam"
        experiment_minus = "{outdir}/bams_stranded/{sample}.{genome}.experiment.m.bam"
        index = ["{outdir}/bams_stranded/{sample}.{genome}.experiment.p.bam.bai",
                 "{outdir}/bams_stranded/{sample}.{genome}.experiment.m.bam.bai"]
        if controls != "None":
            control_plus = "{outdir}/bams_stranded/{sample}.{genome}.control.p.bam"
            control_minus = "{outdir}/bams_stranded/{sample}.{genome}.control.m.bam"
            index.extend(["{outdir}/bams_stranded/{sample}.{genome}.control.p.bam.bai",
                         "{outdir}/bams_stranded/{sample}.{genome}.control.m.bam.bai"]),
            macs2_ss_command="macs2 callpeak -t {input.experiment_plus} -c {input.control_plus} {params.callpeak} --nomodel" \
                             + " --extsize ${{frag_med%.*}} --outdir {wildcards.outdir}/peaks_macs_stranded/ " \
                             + "-n {wildcards.sample}_{wildcards.genome}_plus && " \
                             + "macs2 callpeak -t {input.experiment_minus} -c {input.control_minus} {params.callpeak} --nomodel" \
                             + " --extsize ${{frag_med%.*}} --outdir {wildcards.outdir}/peaks_macs_stranded/ " \
                             + "-n {wildcards.sample}_{wildcards.genome}_minus"
            epic2_ss_command="epic2 -t {input.experiment_plus} -c {input.control_plus} {params.epic} -fs ${{frag_med%.*}}" \
                             + " --output {output.plus} && " \
                             + "epic2 -t {input.experiment_minus} -c {input.control_minus} {params.epic} -fs ${{frag_med%.*}}" \
                             + " --output {output.minus}"
        else:
            control_plus = experiment_plus
            control_minus = experiment_minus
            macs2_ss_command="macs2 callpeak -t {input.experiment_plus} {params.callpeak} --nomodel" \
                             + " --extsize ${{frag_med%.*}} --outdir {wildcards.outdir}/peaks_macs_stranded/ " \
                             + "-n {wildcards.sample}_{wildcards.genome}_plus && " \
                             + "macs2 callpeak -t {input.experiment_minus} {params.callpeak} --nomodel" \
                             + " --extsize ${{frag_med%.*}} --outdir {wildcards.outdir}/peaks_macs_stranded/ " \
                             + "-n {wildcards.sample}_{wildcards.genome}_minus"
            epic2_ss_command="epic2 -t {input.experiment_plus} {params.epic} -fs ${{frag_med%.*}}" \
                             + " --output {output.plus} && " \
                             + "epic2 -t {input.experiment_minus} {params.epic} -fs ${{frag_med%.*}}" \
                             + " --output {output.minus}"

        rule deeptools_coverage_stranded:
                input:
                    bam="{outdir}/bams/{sample}.{genome}.experiment.bam",
                    index="{outdir}/bams/{sample}.{genome}.experiment.bam.bai"
                output:
                    plus="{outdir}/coverage_stranded/{sample}.{genome}.p.bw",
                    minus="{outdir}/coverage_stranded/{sample}.{genome}.m.bw"
                threads: cores
                log: "{outdir}/logs/{sample}_{genome}__deeptools_coverage_stranded.log"
                params:
                    extra="--ignoreForNormalization chrX chrY chrM --ignoreDuplicates --minMappingQuality" \
                          + " 30 --binSize 10 --effectiveGenomeSize " + str(effective_genome_size)
                shell: """
                    (bamCoverage -b {input.bam} -p {threads} {params.extra} --filterRNAstrand forward -o {output.plus}
                    bamCoverage -b {input.bam} -p {threads} {params.extra} --filterRNAstrand reverse -o {output.minus}) &> {log}
                """

        if paired_end:
            rule deeptools_get_pe_fragment_sizes:
                # TODO: First, remove the low MAPQ score reads as they increase the overall size
                input:
                    bam="{outdir}/bams/{sample}.{genome}.experiment.bam",
                    index="{outdir}/bams/{sample}.{genome}.experiment.bam.bai"
                output:
                    info="{outdir}/info/{sample}.{genome}.experiment.frag_lengths.txt",
                    plot="{outdir}/info/{sample}.{genome}.experiment.frag_lengths.png"
                params: "--samplesLabel " + sample_name
                log: "{outdir}/logs/{sample}_{genome}__deeptools_get_pe_fragment_sizes.log"
                threads: cores
                shell: "(bamPEFragmentSize --bamfiles {input.bam} --histogram {output.plot} --table {output.info} -p {threads}) &> {log}"

            rule epic_callpeaks_pe_stranded:
                input:
                    info="{outdir}/info/{sample}.{genome}.experiment.frag_lengths.txt",
                    experiment_plus=experiment_plus,
                    experiment_minus=experiment_minus,
                    control_plus=control_plus,
                    control_minus=control_minus
                output:
                    plus="{outdir}/peaks_epic_stranded/{sample}_{genome}_plus.bed",
                    minus="{outdir}/peaks_epic_stranded/{sample}_{genome}_minus.bed"
                params:
                    epic="-e 100 -fdr .01 --effective-genome-fraction " + str(effective_genome_fraction) \
                         + " --mapq 30 --chromsizes {outdir}/tmp/" + genome + ".chrom.sizes",
                    wget_in="ftp://hgdownload.soe.ucsc.edu/goldenPath/" + genome +\
                         "/bigZips/" + genome + ".chrom.sizes",
                    wget_out="{outdir}/tmp/" + genome + ".chrom.sizes"
                log: "{outdir}/logs/{sample}_{genome}__epic_callpeaks_pe_stranded.log"
                shell: "(frag_med=$(head -n 2 {input.info} | tail -n 1 | awk '{{print $6}}')" \
                       + " && wget {params.wget_in} -O {params.wget_out} && " + epic2_ss_command + ") &> {log}"

            rule macs2_callpeaks_pe_stranded:
                input:
                    info="{outdir}/info/{sample}.{genome}.experiment.frag_lengths.txt",
                    experiment_plus=experiment_plus,
                    experiment_minus=experiment_minus,
                    control_plus=control_plus,
                    control_minus=control_minus
                output:
                    plus="{outdir}/peaks_macs_stranded/{sample}_{genome}_plus_peaks.xls",
                    minus="{outdir}/peaks_macs_stranded/{sample}_{genome}_minus_peaks.xls"
                params:
                    callpeak="--broad-cutoff .01 -q .01 --broad -g " + str(effective_genome_size)
                log: "{outdir}/logs/{sample}_{genome}__macs2_callpeaks_pe_stranded.log"
                shell: "(frag_med=$(head -n 2 {input.info} | tail -n 1 | awk '{{print $6}}') && " \
                        + macs2_ss_command + ") &> {log}"

        else:
            rule epic_callpeaks_se_stranded:
                input:
                    info="{outdir}/info/{sample}.{genome}.experiment_predictd.txt",
                    experiment_plus=experiment_plus,
                    experiment_minus=experiment_minus,
                    control_plus=control_plus,
                    control_minus=control_minus
                output:
                    plus="{outdir}/peaks_epic_stranded/{sample}_{genome}_plus.bed",
                    minus="{outdir}/peaks_epic_stranded/{sample}_{genome}_minus.bed"
                params:
                    epic="-e 100 -fdr .01 --effective-genome-fraction " + str(effective_genome_fraction) \
                         + " --mapq 30 --chromsizes {outdir}/tmp/" + genome + ".chrom.sizes",
                    wget_in="ftp://hgdownload.soe.ucsc.edu/goldenPath/" + genome +\
                         "/bigZips/" + genome + ".chrom.sizes",
                    wget_out="{outdir}/tmp/" + genome + ".chrom.sizes"
                log: "{outdir}/logs/{sample}_{genome}__epic_callpeaks_se_stranded.log"
                shell: "(frag_med=$(cat {input.info} | grep -o 'predicted fragment length is [0-9]* bps' | cut -d ' ' -f 5)" \
                       + " && wget {params.wget_in} -O {params.wget_out} && " + epic2_ss_command + ") &> {log}"

            rule macs2_callpeaks_se_stranded:
                # Based on https://github.com/PEHGP/drippipline
                input:
                    info="{outdir}/info/{sample}.{genome}.experiment_predictd.txt",
                    experiment_plus=experiment_plus,
                    experiment_minus=experiment_minus,
                    control_plus=control_plus,
                    control_minus=control_minus
                output:
                    plus="{outdir}/peaks_macs_stranded/{sample}_{genome}_plus_peaks.xls",
                    minus="{outdir}/peaks_macs_stranded/{sample}_{genome}_minus_peaks.xls"
                params:
                    callpeak="--broad-cutoff .01 -q .01 --broad -g " + str(effective_genome_size)
                log: "{outdir}/logs/{sample}_{genome}__macs2_callpeaks_se_stranded.log"
                shell: "(frag_med=$(cat {input.info} | grep -o 'predicted fragment length is [0-9]* bps' | cut -d ' ' -f 5)" \
                       + " && " + macs2_ss_command + ") &> {log}"


## TODO: Correlation module goes here
rule get_correlation_bin_scores:
    input: "{outdir}/coverage_unstranded/{sample}.{genome}.bw",
    params:
        helpers_dir = helpers_dir
    output:
        npz=temp("{outdir}/QC/{sample}.{genome}.gold_standard_bin_scores.npz"),
        tab="{outdir}/QC/{sample}.{genome}.gold_standard_bin_scores.tab"
    threads: cores
    log: "{outdir}/logs/{sample}_{genome}__get_correlation_bin_scores.log"
    shell:
        """
        (multiBigwigSummary BED-file --BED {params.helpers_dir}/data/correlation_genes_100kb.{wildcards.genome}.1kbwindow.bed \
        -o {output.npz} -b {input} --outRawCounts {output.tab} -p {threads}) &> {log}
        """


rule correlation_analysis:
    input: "{outdir}/QC/{sample}.{genome}.gold_standard_bin_scores.tab"
    output:
        image="{outdir}/QC/{sample}.{genome}.correlation.png",
        data="{outdir}/QC/{sample}.{genome}.correlation.rda"
    params:
        helpers_dir = helpers_dir,
        mode = mode
    log: "{outdir}/logs/{sample}_{genome}__correlation_analysis.log"
    shell:
         """
         (Rscript {params.helpers_dir}/correlation_test.R {input} {wildcards.sample} {params.mode} {params.helpers_dir}) &> {log}
         """


rule assign_genome_annotations:
    # Calculates the percentage of called peaks overlapping with repeat regions
    input:
        gene_anno=genome_home_dir + "/{genome}/homer_anno.txt",
        peaks_macs2="{outdir}/peaks_macs_unstranded/{sample}_{genome}_peaks.broadPeak"
    output:
        stats_out="{outdir}/QC/{sample}_{genome}.feature_overlaps.txt"
    log: "{outdir}/logs/{sample}_{genome}__assign_genome_annotations.log"
    shell: """
        (assignGenomeAnnotation {input.peaks_macs2} {input.gene_anno} > {output.stats_out}) &> {log}
    """


rule prepare_report:
    input: final_report_input
    output:
        html="{outdir}/{sample}_{genome}.QC_report.html",
        rda="{outdir}/{sample}_{genome}.QC_report.rda"
    params:
        helpers_dir=helpers_dir,
        sample_name=sample_name,
        configs=outdirorig + "rseqVars.json",
        final_report_dict_file=final_report_dict_file[0]
    log: "{outdir}/logs/{sample}_{genome}__prepare_report.log"
    shell:
     """
     (echo {output}
     Rscript {params.helpers_dir}/prepare_report.R {params.final_report_dict_file} {params.sample_name} {params.configs}
     rm {params.final_report_dict_file}) &> {log}
     """


rule run_QmRLFS_finder:
    # Packaged with RSeq from: https://github.com/piroonj/QmRLFS-finder (Oct 12 2020)
    input: genome_home_dir + "/{genome}/{genome}.fa"
    output:
        table=genome_home_dir + "/{genome}/rloop_predictions/RLFS.{genome}.out.table.txt",
        bed=genome_home_dir + "/{genome}/rloop_predictions/RLFS.{genome}.out.table.bed"
    params:
        helpers_dir=helpers_dir,
        genome_home_dir=genome_home_dir
    shell: """
    python {params.helpers_dir}/external/QmRLFS-finder.py -i {input} \
    -o {params.genome_home_dir}/{wildcards.genome}/rloop_predictions/RLFS.{wildcards.genome}
    awk '{{OFS="\t";print($1,$4,$14,$3,0,$21)}}' {output.table} > {output.bed}'
    """


rule rlfs_enrichment:
    input:
         peaks="{outdir}/peaks_macs_unstranded/{sample}_{genome}_peaks.broadPeak",
         rlfs=genome_home_dir + "/{genome}/rloop_predictions/RLFS.{genome}.out.table.bed"
    output: "{outdir}/QC/{sample}_{genome}.rlfs_enrichment.rda"
    params:
        helpers_dir=helpers_dir
    threads: cores
    log: "{outdir}/logs/{sample}_{genome}__rlfs_enrichment.log"
    shell: """
     (Rscript {params.helpers_dir}/rlfs_perm_test.R {threads} {wildcards.genome} {input.peaks} {input.rlfs} {output}) &> {log}
     """








