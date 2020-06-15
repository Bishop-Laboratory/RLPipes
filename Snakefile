import json
import os
output_json = json.load(open('run_vars.json'))
# SAMPLES = output_json['SRX1798990_24hr_input_S_96_DRIP_Seq'] # For merging SE reps
# SAMPLES = output_json['SRX2481503_S96_DRIP_Seq_1'] # For paired_end unstranded
SAMPLES = output_json['SRX6427720_DRB_qDRIP-seq_1'] # For paired_end and strand-specific
# SAMPLES = output_json['SRX1070678_NT2_DRIP-seq_1']

# Basic vars
paired_end=SAMPLES['paired_end'][0]
strand_specific=SAMPLES['strand_specific'][0]
genome=SAMPLES['genome'][0]
genome_home_dir=SAMPLES['genome_home_dir'][0]
effective_genome_size=SAMPLES['effective_genome_size'][0]
full_genome_length=SAMPLES['full_genome_length'][0]
effective_genome_fraction=effective_genome_size/full_genome_length
cores=SAMPLES['cores'][0]
sample_name=SAMPLES['sample_name'][0]
outdir=SAMPLES['out_dir'][0]
experiments=SAMPLES['experiments'][0]
controls=SAMPLES['controls'][0]

if paired_end:
    output_fastq_1 = expand("{outdir}/fastqs_cleaned/{sample_name}_R1.fastq",
                     sample_name=sample_name, outdir=outdir)
    output_fastq_2 = expand("{outdir}/fastqs_cleaned/{sample_name}_R2.fastq",
                     sample_name=sample_name, outdir=outdir)
    merge_input_1 = expand("{outdir}/fastqs_raw/{sample_name}/{srr_acc}.sra_1.fastq",
                           sample_name=sample_name,
                           srr_acc=experiments, outdir=outdir)
    merge_input_2 = expand("{outdir}/fastqs_raw/{sample_name}/{srr_acc}.sra_2.fastq",
                           sample_name=sample_name,
                           srr_acc=experiments, outdir=outdir)
else:
    output_fastq = expand("{outdir}/fastqs_cleaned/{sample_name}.fastq",
                     sample_name=sample_name, outdir=outdir)
    merge_input = expand("{outdir}/fastqs_raw/{sample_name}/{srr_acc}.sra.fastq",
                         sample_name=sample_name,
                         srr_acc=experiments, outdir=outdir)

if strand_specific:
    coverage_output = expand("{outdir}/coverage/{sample}.{genome}.p.bw", genome=genome,
                        sample=sample_name, outdir=outdir)
else:
    coverage_output = expand("{outdir}/coverage/{sample}.{genome}.bw", genome=genome,
                            sample=sample_name, outdir=outdir)

rule output:
        input:
            # expand("{outdir}/libtype_infer/{sample_name}.{genome}/lib_format_counts.json", genome=SAMPLES['genome'],
            #        sample_name=SAMPLES['sample_name'], outdir=SAMPLES['out_dir'])
            #  expand("{outdir}/bams_stranded/{sample}.{genome}.p.bam", genome=SAMPLES['genome'],
            #        sample=SAMPLES['sample_name'], outdir=SAMPLES['out_dir']),
            #  expand("{outdir}/coverage/{sample}.{genome}.p.bw", genome=SAMPLES['genome'],
            #        sample=SAMPLES['sample_name'], outdir=SAMPLES['out_dir']),
             expand("{outdir}/peaks_macs/{sample}_{genome}_peaks.broadPeak", genome=genome,
                   sample=sample_name, outdir=outdir),
             expand("{outdir}/peaks_epic/{sample}_{genome}.bed", genome=genome,
                   sample=sample_name, outdir=outdir),
             expand("{outdir}/peaks_macs/{sample}_{genome}_plus_peaks.xls", genome=genome,
                   sample=sample_name, outdir=outdir),
             coverage_output,
             expand("{outdir}/bams/{sample}.{genome}.bam.bai", genome=genome,
                   sample=sample_name, outdir=outdir),
            # expand(SAMPLES['genome_dir'][0] + "/{genome}/salmon_index/info.json", genome=SAMPLES['genome'],
            #       sample_name=SAMPLES['sample_name'], outdir=SAMPLES['out_dir'])


rule download_sra:
    output: "{outdir}/fastqs_raw/{sample_name}/{srr_acc}.sra"
    shell:
        "prefetch {wildcards.srr_acc} --output-file {output}"

if paired_end:
    rule sra_to_fastq:
        output:
            "{outdir}/fastqs_raw/{sample_name}/{srr_acc}.sra_1.fastq",
            "{outdir}/fastqs_raw/{sample_name}/{srr_acc}.sra_2.fastq"
        input: "{outdir}/fastqs_raw/{sample_name}/{srr_acc}.sra"
        threads: SAMPLES['cores'][0]
        shell:
            "fasterq-dump -e {threads} --split-files -O {wildcards.outdir}/fastqs_raw/{wildcards.sample_name}/ {input}"

    rule merge_fastq:
        output:
            R1="{outdir}/fastqs/{sample_name}_R1.fastq",
            R2="{outdir}/fastqs/{sample_name}_R2.fastq"
        input:
            R1=merge_input_1,
            R2=merge_input_2
        shell: """
            cat {input.R1} > {output.R1}
            cat {input.R2} > {output.R2}
        """

    rule fastp_pe:
        input:
            sample=["{outdir}/fastqs/{sample}_R1.fastq", "{outdir}/fastqs/{sample}_R2.fastq"]
        output:
            trimmed=["{outdir}/fastqs_cleaned/{sample}_R1.fastq", "{outdir}/fastqs_cleaned/{sample}_R2.fastq"],
            html="{outdir}/QC/fastq/{sample}.html",
            json="{outdir}/QC/fastq/{sample}.json"
        log:
            "{outdir}/logs/{sample}_fastp.log"
        params:
            extra=""
        threads: cores
        wrapper:
            "0.60.0/bio/fastp"
else:
    rule sra_to_fastq:
        output: "{outdir}/fastqs_raw/{sample_name}/{srr_acc}.sra.fastq"
        input: "{outdir}/fastqs_raw/{sample_name}/{srr_acc}.sra"
        threads: cores
        shell:
            "fasterq-dump -e {threads} --split-files -O {wildcards.outdir}/fastqs_raw/{wildcards.sample_name}/ {input}"

    rule merge_fastq:
        output: "{outdir}/fastqs/{sample_name}.fastq"
        input: merge_input
        shell: "cat {input} > {output}"

    rule fastp_se:
        input:
            sample=["{outdir}/fastqs/{sample}.fastq"]
        output:
            trimmed="{outdir}/fastqs_cleaned/{sample}.fastq",
            html="{outdir}/QC/fastq/{sample}.html",
            json="{outdir}/QC/fastq/{sample}.json"
        log:
            "{outdir}/logs/{sample}_fastp.log"
        params:
            extra=""
        threads: cores
        wrapper:
            "0.60.0/bio/fastp"

if paired_end:
    bwa_reads=["{outdir}/fastqs_cleaned/{sample}_R1.fastq", "{outdir}/fastqs_cleaned/{sample}_R2.fastq"]
else:
    bwa_reads=["{outdir}/fastqs_cleaned/{sample}.fastq"]

rule download_fasta:
    output:
        genome_home_dir + "/{genome}/{genome}.fa"
    shell: """
        wget -O {output}.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz
        gunzip {output}.gz
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
        "0.60.0/bio/bwa/index"

rule bwa_mem:
    input:
        bwa_index_done=genome_home_dir + "/{genome}/bwa_index/{genome}.amb",
        reads=bwa_reads
    output:
        bam="{outdir}/bam_dups/{sample}.{genome}.bam"
    log:
        "{outdir}/logs/{sample}_{genome}_bwa.log"
    params:
        index=genome_home_dir + "/{genome}/bwa_index/{genome}",
        bwa_extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        picard_extra="",
        samtools_sort_extra=""
    threads: cores
    shell: """
        bwa mem -t {threads} {params.bwa_extra} {params.index} {input.reads} | \
        samtools view -b -@ {threads} - | \
        samtools sort {params.samtools_sort_extra} -@ {threads} -o {output.bam} - 
     """

rule mark_duplicates:
    input: "{outdir}/bam_dups/{sample}.{genome}.bam"
    output:
        bam="{outdir}/bams/{sample}.{genome}.bam",
        metrics="{outdir}/info/{sample}.{genome}.duplicate.metrics.txt"
    params:
        ""
    threads: cores
    shell: """
        picard MarkDuplicates I={input} O={output.bam} M={output.metrics}
    """

rule index_bam:
    input: "{outdir}/bams/{sample}.{genome}.bam"
    output: "{outdir}/bams/{sample}.{genome}.bam.bai"
    threads: cores
    shell: """
        samtools index -@ {threads} {input}
    """

if strand_specific and paired_end:
    rule split_strands_pe:
        input: "{outdir}/bams/{sample}.{genome}.bam",
        output:
            plus="{outdir}/bams_stranded/{sample}.{genome}.p.bam",
            minus="{outdir}/bams_stranded/{sample}.{genome}.m.bam"
        threads: cores
        shell: """
            # Adapted from https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
            (mkdir {wildcards.outdir}/bams_tmp) &> /dev/null || true
            samtools view -@ {threads} -b -f 128 -F 16 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.fwd1.bam
            samtools view -@ {threads} -b -f 80 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.fwd2.bam
            samtools merge -@ {threads} -f {output.plus} {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.fwd1.bam {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.fwd2.bam
            samtools index -@ {threads} {output.plus}
            samtools view -@ {threads} -b -f 144 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.rev1.bam
            samtools view -@ {threads} -b -f 64 -F 16 {input} > {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.rev2.bam
            samtools merge -@ {threads} -f {output.minus} {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.rev1.bam {wildcards.outdir}/bams_tmp/{wildcards.sample}.{wildcards.genome}.rev2.bam
            samtools index -@ {threads} {output.minus}
        """
elif strand_specific:
    rule split_strands_se:
        input: "{outdir}/bams/{sample}.{genome}.bam",
        output:
            plus="{outdir}/bams_stranded/{sample}.{genome}.p.bam",
            minus="{outdir}/bams_stranded/{sample}.{genome}.m.bam"
        threads: cores
        shell: """
            # Adapted from https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
            (mkdir {wildcards.outdir}/bams_tmp) &> /dev/null || true
            samtools view -@ {threads} -b -F 16 {input} > {output.plus}
            samtools view -@ {threads} -b -f 16 {input} > {output.minus}
        """

if paired_end:
    macs2_params="--broad -g -f BAMPE " + str(effective_genome_size)
else:
    macs2_params="--broad -g " + str(effective_genome_size)

rule macs2_callpeak:
    input:
        "{outdir}/bams/{sample}.{genome}.bam",
    output:
        "{outdir}/peaks_macs/{sample}_{genome}_peaks.broadPeak",
    params: macs2_params
    shell: """
        macs2 callpeak -t {input} {params} --outdir {wildcards.outdir}/peaks_macs/ -n {wildcards.sample}_{wildcards.genome}
    """

if paired_end:
    rule bampe_to_bedpe:
        input:
            "{outdir}/bams/{sample}.{genome}.bam"
        output:
            "{outdir}/bedpes/{sample}.{genome}.bed"
        shell:


    rule epic2_pe_callpeak:
        input:
            "{outdir}/bedpes/{sample}.{genome}.bed",
        output:
            "{outdir}/peaks_epic/{sample}_{genome}.bed",
        params:
            epic="--effective-genome-fraction " + str(effective_genome_fraction) \
                 + " --chromsizes {outdir}/tmp/" + genome + ".chrom.sizes --output peaks.txt",
            wget_in="ftp://hgdownload.soe.ucsc.edu/goldenPath/" + genome +\
                 "/bigZips/" + SAMPLES['genome'][0] + ".chrom.sizes",
            wget_out="{outdir}/tmp/" + genome + ".chrom.sizes"
        shell: """
            wget {params.wget_in} -O {params.wget_out}
            epic2 --treatment {input} {params.epic} --output {output}
        """
else:
    rule epic2_se_callpeak:
        input:
            "{outdir}/bams/{sample}.{genome}.bam",
        output:
            "{outdir}/peaks_epic/{sample}_{genome}.bed",
        params:
            epic="--effective-genome-fraction " + str(effective_genome_fraction) \
                 + " --chromsizes {outdir}/tmp/" + genome + ".chrom.sizes --output peaks.txt",
            wget_in="ftp://hgdownload.soe.ucsc.edu/goldenPath/" + genome +\
                 "/bigZips/" + SAMPLES['genome'][0] + ".chrom.sizes",
            wget_out="{outdir}/tmp/" + genome + ".chrom.sizes"
        shell: """
            wget {params.wget_in} -O {params.wget_out}
            epic2 --treatment {input} {params.epic} --output {output}
        """



if strand_specific and paired_end:
    rule deeptools_get_pe_fragment_sizes:
        input: "{outdir}/bams/{sample}.{genome}.bam",
        output:
            info="{outdir}/info/{sample}.{genome}.frag_lengths.txt",
            plot="{outdir}/info/{sample}.{genome}.frag_lengths.png"
        params: "--samplesLabel " + sample_name
        threads: cores
        shell: "bamPEFragmentSize --bamfiles {input} --histogram {output.plot} --table {output.info} -p {threads}"


    # rule epic_callpeaks_pe_stranded:
    #     input:
    #         info="{outdir}/info/{sample}.{genome}.frag_lengths.txt",
    #         plus="{outdir}/bams_stranded/{sample}.{genome}.p.bam",
    #         minus="{outdir}/bams_stranded/{sample}.{genome}.m.bam"
    #     output:
    #         plus="{outdir}/peaks_macs/{sample}_{genome}_plus_peaks.xls",
    #         minus="{outdir}/peaks_macs/{sample}_{genome}_minus_peaks.xls"
    #     params:
    #         epic="--effective-genome-fraction " + str(effective_genome_fraction) \
    #              + " --chromsizes {outdir}/tmp/" + genome + ".chrom.sizes --output peaks.txt",
    #         wget_in="ftp://hgdownload.soe.ucsc.edu/goldenPath/" + genome +\
    #              "/bigZips/" + SAMPLES['genome'][0] + ".chrom.sizes",
    #         wget_out="{outdir}/tmp/" + genome + ".chrom.sizes"
    #     shell: """
    #         frag_med=$(head -n 2 {input.info} | tail -n 1 | awk "{{print \$6}}")
    #         frag_min=$(head -n 2 {input.info} | tail -n 1 | awk "{{print \$11}}")
    #         frag_max=$(head -n 2 {input.info} | tail -n 1 | awk "{{print \$19}}")
    #         wget {params.wget_in} -O {params.wget_out}
    #         epic2 --treatment {input} {params.epic} --output {output}
    #
    #         macs2 callpeak -t {input.plus} {params.extra} --min-length $frag_min --nomodel --ext-size $frag_med --outdir {wildcards.outdir}/peaks_macs/ -n {wildcards.sample}_{wildcards.genome}_plus
    #         macs2 callpeak -t {input.minus} {params.extra} --min-length $frag_min --nomodel --ext-size $frag_med --outdir {wildcards.outdir}/peaks_macs/ -n {wildcards.sample}_{wildcards.genome}_minus
    #     """

    rule macs2_callpeaks_pe_stranded:
        input:
            info="{outdir}/info/{sample}.{genome}.frag_lengths.txt",
            plus="{outdir}/bams_stranded/{sample}.{genome}.p.bam",
            minus="{outdir}/bams_stranded/{sample}.{genome}.m.bam"
        output:
            plus="{outdir}/peaks_macs/{sample}_{genome}_plus_peaks.xls",
            minus="{outdir}/peaks_macs/{sample}_{genome}_minus_peaks.xls"
        params:
            callpeak="--broad -g " + str(effective_genome_size)
        shell: """
            frag_med=$(head -n 2 {input.info} | tail -n 1 | awk "{{print \$6}}")
            frag_min=$(head -n 2 {input.info} | tail -n 1 | awk "{{print \$11}}")
            frag_max=$(head -n 2 {input.info} | tail -n 1 | awk "{{print \$19}}")
            macs2 callpeak -t {input.plus} {params.callpeak} --min-length $frag_min --nomodel --ext-size $frag_med --outdir {wildcards.outdir}/peaks_macs/ -n {wildcards.sample}_{wildcards.genome}_plus
            macs2 callpeak -t {input.minus} {params.callpeak} --min-length $frag_min --nomodel --ext-size $frag_med --outdir {wildcards.outdir}/peaks_macs/ -n {wildcards.sample}_{wildcards.genome}_minus
        """

    rule deeptools_coverage_stranded:
        input:
            bam="{outdir}/bams/{sample}.{genome}.bam",
            index="{outdir}/bams/{sample}.{genome}.bam.bai"
        output:
            plus="{outdir}/coverage/{sample}.{genome}.p.bw",
            minus="{outdir}/coverage/{sample}.{genome}.m.bw"
        threads: cores
        params:
            extra="--ignoreForNormalization chrX chrY chrM --ignoreDuplicates --minMappingQuality" \
                  + " 30 --binSize 20 --effectiveGenomeSize " + str(effective_genome_size)
        shell: """
            bamCoverage -b {input.bam} -p {threads} {params.extra} --filterRNAstrand forward -o {output.plus}
            bamCoverage -b {input.bam} -p {threads} {params.extra} --filterRNAstrand reverse -o {output.minus}
        """
elif strand_specific:
    # TODO: Need to write single-end stranded caller
    rule macs2_predictd:
        input: "{outdir}/bams/{sample}.{genome}.bam"
        output: "{outdir}/info/{sample}.{genome}.macs2_d.txt"
        params: "-g " + str(effective_genome_size) + " --outdir " + outdir + "/info"
        shell: "macs2 predictd -i {input} > {output}"

    rule macs2_callpeaks_se_stranded:
        # Based on https://github.com/PEHGP/drippipline
        input:
            info="{outdir}/info/{sample}.{genome}.macs2_d.txt",
            plus="{outdir}/bams_stranded/{sample}.{genome}.p.bam",
            minus="{outdir}/bams_stranded/{sample}.{genome}.m.bam"
        output:
            plus="{outdir}/peaks/{sample}_{genome}_plus_peaks.xls",
            minus="{outdir}/peaks/{sample}_{genome}_minus_peaks.xls"
        params:
            predictd="-f BAMPE -g " + str(effective_genome_size),
            callpeak="--broad -f BAMPE -g " + str(effective_genome_size)
        shell: """
            frag_size=$(cat {input.info} | grep 'predicted fragment length is ([0-9]+) bps')
            macs2 callpeak -t {input.plus} {params.callpeak} --no-model --ext-size $frag_size --outdir {wildcards.outdir}/peaks/ -n {wildcards.sample}_{wildcards.genome}_plus
            macs2 callpeak -t {input.minus} {params.callpeak} --no-model --ext-size $frag_size --outdir {wildcards.outdir}/peaks/ -n {wildcards.sample}_{wildcards.genome}_minus
        """
else:
    rule deeptools_coverage_unstranded:
        input:
            bam="{outdir}/bams/{sample}.{genome}.bam",
            index="{outdir}/bams/{sample}.{genome}.bam.bai"
        output: "{outdir}/coverage/{sample}.{genome}.bw"
        threads: cores
        params:
            extra="--ignoreForNormalization chrX chrY chrM --ignoreDuplicates --minMappingQuality" \
                  + " 30 --binSize 20 --effectiveGenomeSize " + str(effective_genome_size)
        shell: """
            bamCoverage -b {input.bam} -p {threads} {params.extra} -o {output}
        """



rule download_annotation:
    output:
        genome_home_dir + "/{genome}/{genome}.gtf"
    shell: """
        wget -O {output}.gz \
        ftp://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/genes/{wildcards.genome}.refGene.gtf.gz
        gunzip {output}.gz
    """


