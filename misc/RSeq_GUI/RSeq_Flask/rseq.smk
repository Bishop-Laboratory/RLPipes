########################################################################################################################
##############################################   Parse inputs    #######################################################
########################################################################################################################

# Basic vars
paired_end=config['paired_end'][0]
strand_specific=config['strand_specific'][0]
genome=config['genome'][0]
genome_home_dir=config['genome_home_dir'][0]
effective_genome_size=config['effective_genome_size'][0]
full_genome_length=config['full_genome_length'][0]
effective_genome_fraction = effective_genome_size/full_genome_length
cores=config['cores'][0]
sample_name=config['sample_name'][0]
outdir=config['out_dir'][0]
experiments=config['experiments'][0]
controls=config['controls'][0]


# Set expected output and inputs for merging of technical replicates
output_fastq_experiment_1 = expand("{outdir}/fastqs/{sample_name}_experiment_R1.fastq",
                 sample_name=sample_name, outdir=outdir)
output_fastq_experiment_2 = expand("{outdir}/fastqs/{sample_name}_experiment_R2.fastq",
                 sample_name=sample_name, outdir=outdir)
merge_input_experiment_1 = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}.sra_1.fastq",
                           sample_name=sample_name,
                           srr_acc=experiments, outdir=outdir)
merge_input_experiment_2 = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}.sra_2.fastq",
                           sample_name=sample_name,
                           srr_acc=experiments, outdir=outdir),
output_fastq_control_1 = expand("{outdir}/fastqs/{sample_name}_control_R1.fastq",
                 sample_name=sample_name, outdir=outdir)
output_fastq_control_2 = expand("{outdir}/fastqs/{sample_name}_control_R2.fastq",
                 sample_name=sample_name, outdir=outdir)
merge_input_control_1 = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}.sra_1.fastq",
                           sample_name=sample_name,
                           srr_acc=controls, outdir=outdir)
merge_input_control_2 = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}.sra_2.fastq",
                           sample_name=sample_name,
                           srr_acc=controls, outdir=outdir)
output_fastq_experiment = expand("{outdir}/fastqs/{sample_name}.fastq",
                 sample_name=sample_name, outdir=outdir)
merge_input_experiment = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}.sra.fastq",
                         sample_name=sample_name,
                         srr_acc=experiments, outdir=outdir),
output_fastq_control = expand("{outdir}/fastqs/{sample_name}.fastq",
                 sample_name=sample_name, outdir=outdir)
merge_input_control = expand("{outdir}/tmp/fastqs_raw/{sample_name}/{srr_acc}.sra.fastq",
                         sample_name=sample_name,
                         srr_acc=controls, outdir=outdir)

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

# TODO: Set the TMPDIR

########################################################################################################################
##############################################   Main pipeline    ######################################################
########################################################################################################################

rule output:
        input:
            peak_output,
            coverage_output

rule download_sra:
    output: temp("{outdir}/tmp/sras/{sample}/{srr_acc}.sra")
    log: "{outdir}/logs/{sample}/{srr_acc}__download_sra.log"
    shell: "(prefetch {wildcards.srr_acc} --output-file {output}) &> {log}"

if paired_end:
    rule sra_to_fastq_pe:
        input:
            "{outdir}/tmp/sras/{sample}/{srr_acc}.sra"
        output:
            temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.sra_1.fastq"),
            temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.sra_2.fastq")
        threads: cores
        log: "{outdir}/logs/{sample}/{srr_acc}__sra_to_fastq_pe.log"
        shell:
            "(fasterq-dump -e {threads} --split-files -O" \
            + " {wildcards.outdir}/tmp/fastqs_raw/{wildcards.sample}/ {input}) &> {log}"

    rule merge_fastq_pe_experiment:
        input:
            R1=merge_input_experiment_1,
            R2=merge_input_experiment_2
        output:
            R1=temp("{outdir}/tmp/fastqs_dup/{sample}_experiment_R1.fastq"),
            R2=temp("{outdir}/tmp/fastqs_dup/{sample}_experiment_R2.fastq")
        log: "{outdir}/logs/{sample}/merge_fastq_pe_experiment.log"
        shell: """
            ( cat {input.R1} > {output.R1}
            cat {input.R2} > {output.R2} ) &> {log}
        """

    rule merge_fastq_pe_control:
        input:
            R1=merge_input_control_1,
            R2=merge_input_control_2
        output:
            R1=temp("{outdir}/tmp/fastqs_dup/{sample}_control_R1.fastq"),
            R2=temp("{outdir}/tmp/fastqs_dup/{sample}_control_R2.fastq")
        log: "{outdir}/logs/{sample}/merge_fastq_pe_control.log"
        shell: """
            ( cat {input.R1} > {output.R1}
            cat {input.R2} > {output.R2} ) &> {log}
        """

    rule clumpify_pe:
        input:
            R1="{outdir}/tmp/fastqs_dup/{sample}_{exp_type}_R1.fastq",
            R2="{outdir}/tmp/fastqs_dup/{sample}_{exp_type}_R2.fastq"
        output:
            R1=temp("{outdir}/tmp/fastqs_dedup/{sample}_{exp_type}_R1.fastq"),
            R2=temp("{outdir}/tmp/fastqs_dedup/{sample}_{exp_type}_R2.fastq")
        params: ""
        threads: cores
        log: "{outdir}/logs/{sample}/{exp_type}__clumpify_pe.log"
        shell:
            "(clumpify.sh t={threads} {params} in1={input.R1}" \
            + " in2={input.R2} out1={output.R1} out2={output.R2} dedupe) &> {log}"

    rule fastp_pe:
        input:
            sample=["{outdir}/tmp/fastqs_dedup/{sample}_{exp_type}_R1.fastq",
                    "{outdir}/tmp/fastqs_dedup/{sample}_{exp_type}_R2.fastq"]
        output:
            trimmed=["{outdir}/fastqs/{sample}_{exp_type}_R1.fastq",
                     "{outdir}/fastqs/{sample}_{exp_type}_R2.fastq"],
            html="{outdir}/QC/fastq/{sample}.{exp_type}.html",
            json="{outdir}/QC/fastq/{sample}.{exp_type}.json"
        log: "{outdir}/logs/{sample}/{exp_type}__fastp_pe.log"
        params:
            extra=""
        threads: cores
        wrapper:
            "0.63.0/bio/fastp"

else:

    rule sra_to_fastq_se:
        input: "{outdir}/tmp/sras/{sample}/{srr_acc}.sra"
        output: temp("{outdir}/tmp/fastqs_raw/{sample}/{srr_acc}.sra.fastq")
        threads: cores
        log: "{outdir}/logs/{sample}/{srr_acc}__sra_to_fastq_se.log"
        shell:
            "(fasterq-dump -e {threads} --split-files -O" \
            + " {wildcards.outdir}/tmp/fastqs_raw/{wildcards.sample}/ {input}) &> {log}"

    rule merge_fastq_se_experiment:
        input: merge_input_experiment
        output: temp("{outdir}/tmp/fastqs_dup/{sample}_experiment.fastq")
        log: "{outdir}/logs/{sample}/merge_fastq_se_experiment.log"
        shell: "(cat {input} > {output}) &> {log}"

    rule merge_fastq_se_control:
        input: merge_input_control
        output: temp("{outdir}/tmp/fastqs_dup/{sample}_control.fastq")
        log: "{outdir}/logs/{sample}/merge_fastq_se_control.log"
        shell: "(cat {input} > {output}) &> {log}"

    rule clumpify_se:
        input: "{outdir}/tmp/fastqs_dup/{sample}_{exp_type}.fastq"
        output: temp("{outdir}/tmp/fastqs_dedup/{sample}_{exp_type}.fastq")
        params: ""
        threads: cores
        log: "{outdir}/logs/{sample}/{exp_type}__clumpify_se.log"
        shell: "(clumpify.sh t={threads} {params} in={input} out={output} dedupe) &> {log}"

    rule fastp_se:
        input:
             sample=["{outdir}/tmp/fastqs_dedup/{sample}_{exp_type}.fastq"]
        output:
            trimmed="{outdir}/fastqs/{sample}_{exp_type}.fastq",
            html="{outdir}/QC/fastq/{sample}.{exp_type}.html",
            json="{outdir}/QC/fastq/{sample}.{exp_type}.json"
        log: "{outdir}/logs/{sample}/{exp_type}__fastp_se.log"
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
          check=outdir + "/logs/" + sample_name + "/{genome}__download_fasta.log"
    shell: """
        (wget -O {output}.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz
        gunzip {output}.gz && echo "Indexing BAM at {params.prefix} ... This may take some time ...") &> {params.check}
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
    log: "{outdir}/logs/{sample}/{genome}_{exp_type}__bwa_mem.log"
    params:
        index=genome_home_dir + "/{genome}/bwa_index/{genome}",
        bwa_extra=r"-R '@RG\tID:{sample}_{exp_type}\tSM:{sample}_{exp_type}'",
        picard_extra="",
        samtools_sort_extra=""
    threads: cores
    shell: """
        (bwa mem -t {threads} {params.bwa_extra} {params.index} {input.reads} | \
        samtools view -b -@ {threads} - | \
        samtools sort {params.samtools_sort_extra} -@ {threads} -o {output.bam} -) &> {log}
     """

rule index_bam:
    input: "{outdir}/bams/{sample}.{genome}.{exp_type}.bam"
    output: "{outdir}/bams/{sample}.{genome}.{exp_type}.bam.bai"
    threads: cores
    log: "{outdir}/logs/{sample}/{genome}_{exp_type}__index_bam.log"
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
        log: "{outdir}/logs/{sample}/{genome}_{exp_type}__split_strands_pe.log"
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
        log: "{outdir}/logs/{sample}/{genome}_{exp_type}__split_strands_se.log"
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
            log: "{outdir}/logs/{sample}/{genome}_{exp_type}__macs2_predictd.log"
            shell: "(macs2 predictd -i {input} --outdir {wildcards.outdir}/tmp &> {output}) &> {log}"

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

if controls is not None:
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
    log: "{outdir}/logs/{sample}/{genome}__macs2_callpeak_unstranded.log"
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
    log: "{outdir}/logs/{sample}/{genome}_{exp_type}__bampe_to_bedpe.log"
    shell: """
        (samtools sort {params.samtools_sort_extra} -@ {threads} -o /dev/stdout {input} | \
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
        "{outdir}/peaks_epic_unstranded/{sample}_{genome}.bed",
    params:
        epic="-e 100 -fdr .01 --effective-genome-fraction " + str(effective_genome_fraction) \
             + " --mapq 30 --chromsizes {outdir}/tmp/" + genome + ".chrom.sizes",
        wget_in="ftp://hgdownload.soe.ucsc.edu/goldenPath/" + genome +\
             "/bigZips/" + genome + ".chrom.sizes",
        wget_out="{outdir}/tmp/" + genome + ".chrom.sizes"
    log: "{outdir}/logs/{sample}/{genome}__epic2_callpeak_unstranded.log"
    shell: epic2_command_lt

# Calculate this for all reads
rule deeptools_coverage_unstranded:
        input:
            bam="{outdir}/bams/{sample}.{genome}.experiment.bam",
            index="{outdir}/bams/{sample}.{genome}.experiment.bam.bai"
        output: "{outdir}/coverage_unstranded/{sample}.{genome}.bw"
        threads: cores
        log: "{outdir}/logs/{sample}/{genome}__deeptools_coverage_unstranded.log"
        params:
            extra="--ignoreForNormalization chrX chrY chrM --ignoreDuplicates --minMappingQuality" \
                  + " 30 --binSize 20 --effectiveGenomeSize " + str(effective_genome_size)
        shell: """
            (bamCoverage -b {input.bam} -p {threads} {params.extra} -o {output}) &> {log}
        """

# Handles stranded peak calls and coverage
if strand_specific:
    experiment_plus = "{outdir}/bams_stranded/{sample}.{genome}.experiment.p.bam"
    experiment_minus = "{outdir}/bams_stranded/{sample}.{genome}.experiment.m.bam"
    index = ["{outdir}/bams_stranded/{sample}.{genome}.experiment.p.bam.bai",
             "{outdir}/bams_stranded/{sample}.{genome}.experiment.m.bam.bai"]
    if controls is not None:
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
            log: "{outdir}/logs/{sample}/{genome}__deeptools_coverage_stranded.log"
            params:
                extra="--ignoreForNormalization chrX chrY chrM --ignoreDuplicates --minMappingQuality" \
                      + " 30 --binSize 20 --effectiveGenomeSize " + str(effective_genome_size)
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
            log: "{outdir}/logs/{sample}/{genome}__deeptools_get_pe_fragment_sizes.log"
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
            log: "{outdir}/logs/{sample}/{genome}__epic_callpeaks_pe_stranded.log"
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
            log: "{outdir}/logs/{sample}/{genome}__macs2_callpeaks_pe_stranded.log"
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
            log: "{outdir}/logs/{sample}/{genome}__epic_callpeaks_se_stranded.log"
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
            log: "{outdir}/logs/{sample}/{genome}__macs2_callpeaks_se_stranded.log"
            shell: "(frag_med=$(cat {input.info} | grep -o 'predicted fragment length is [0-9]* bps' | cut -d ' ' -f 5)" \
                   + " && " + macs2_ss_command + ") &> {log}"


rule download_annotation:
    output:
        genome_home_dir + "/{genome}/{genome}.gtf"
    shell: """
        wget -O {output}.gz \
        ftp://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/genes/{wildcards.genome}.refGene.gtf.gz
        gunzip {output}.gz
    """


