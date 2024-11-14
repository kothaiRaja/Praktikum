nextflow.enable.dsl = 2


// Define all parameters
params.reads = "$baseDir/data/reads/rep1_{1,2}.fq.gz"
params.genome = "$launchDir/data/genome.fa"
params.variants = "$launchDir/data/subset_fixed.vcf.gz"
params.denylist = "$launchDir/data/denylist.bed"
params.accession = ""
params.out = "$baseDir/output"

// Output Directory
params.outdir = "$baseDir/output"



process prefetch {
    container "https://depot.galaxyproject.org/singularity/sra-tools%3A2.11.0--pl5321ha49a11a_3"
    publishDir params.cache_dir, mode: "copy"

    input:
    val accession

    output:
    path "${accession}"

    script:
    """
    prefetch ${accession}
    """
}


process fasterqdump {
    container "https://depot.galaxyproject.org/singularity/sra-tools%3A2.11.0--pl5321ha49a11a_3"
    publishDir params.cache_dir, mode: "copy"

    input:
    path sradir

    output:
    path "${sradir}_*.fastq"

    script:
    """
    fasterq-dump --split-files ${sradir}
    """
}


process fastqc_raw {
    tag "${sample_id}"
	container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    publishDir "${params.outdir}/fastqc/raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.zip"), emit: fastqc_zip
    path("*.html"), emit: fastqc_html

    script:
    """
    fastqc --outdir . ${reads.join(' ')}
    """
}
// Fastp Process
process fastp {
    tag "${sample_id}"
	container "https://depot.galaxyproject.org/singularity/fastp%3A0.23.4--hadf994f_3"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.html"), path("${sample_id}_fastp.json"), emit: fastp_reports

    script:
    """
    fastp \
        -i ${reads[0]} -I ${reads[1]} \
        -o ${sample_id}_R1_trimmed.fastq.gz \
        -O ${sample_id}_R2_trimmed.fastq.gz \
		--thread ${task.cpus} \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json
    """
}
process fastqc_trimmed {
    tag "${sample_id}"
    publishDir "${params.outdir}/fastqc/trimmed", mode: 'copy'
	container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.zip"), path("*.html")

    script:
    """
    fastqc --outdir . ${reads.join(' ')}
    """
}



//Create a fasta Genome Index

process PREPARE_GENOME_SAMTOOLS {
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "${params.outdir}/genome/index", mode: "copy"

    input:
    path genome

    output:
    path "${genome}.fai"

    script:
    """
    samtools faidx ${genome}
    """
}

//Create genome dictionary with PICARD for GATK4

process PREPARE_GENOME_PICARD {
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/genome/dict", mode: "copy"

    input:
    path genome

    output:
    path "${genome.baseName}.dict"

    script:
    """
    gatk CreateSequenceDictionary -R $genome -O ${genome.baseName}.dict
   
	# Correct any UR:file path in the resulting genome.dict file to absolute path
    sed -i 's|UR:file:.*|UR:file:/home/kothai/cq-git-sample/Praktikum/data/genome.fa|' ${genome.baseName}.dict
	"""
}

//Create a genome Index file

process PREPARE_STAR_GENOME_INDEX {
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.10b--h6b7c446_1"
    publishDir "${params.outdir}/genome/starindex", mode: "copy"

    input:
    path genome

    output:
    path "genome_index" 

    script:
    """
    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         --genomeDir genome_index \
         --genomeFastaFiles ${genome} \
         --genomeSAindexNbases 11
    """
}



//Filtered and recorded variants

process PREPARE_VCF_FILE {
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.outdir}/genome/vcf", mode: "copy"
	 
    input: 
        path variantsFile
        path denylisted

    output:
        tuple path("${variantsFile.baseName}.filtered.recode.vcf.gz"),
              path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")
  
    script:  
    """
    # Exclude regions from the blacklist using bcftools
    bcftools view -T ^${denylisted} ${variantsFile} -Oz -o ${variantsFile.baseName}.filtered.recode.vcf.gz

    # Index the filtered VCF
    tabix -p vcf ${variantsFile.baseName}.filtered.recode.vcf.gz
    """
}



process RNASEQ_MAPPING_STAR {
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
    publishDir "${params.outdir}/star", mode: "copy"

    input:
    path genomeDir
    tuple val(replicateId), path(reads)

    output:
    tuple val(replicateId), path("Aligned.*.sortedByCoord.out.bam")


    script:
    """
    # Align reads to genome
    STAR --genomeDir $genomeDir \
         --readFilesIn ${reads[0]} ${reads[1]} \
         --runThreadN ${task.cpus} \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outSAMtype BAM SortedByCoordinate \
		 --outSAMattrRGline ID:$replicateId LB:library PL:illumina PU:machine SM:GM12878 \
         --outFileNamePrefix Aligned.
    """
}

process SAMTOOLS_FLAGSTAT {
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/flagstat", mode: "copy"

    input:
    tuple val(replicateId), path(bam)

    output:
    path "${replicateId}_flagstat.txt"

    script:
    """
    samtools flagstat ${bam} > ${replicateId}_flagstat.txt
    """
}


process SAMTOOLS_FILTER_INDEX {
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/star/index", mode: "copy"

    input:
    tuple val(replicateId), path(bam)

    output:
    tuple val(replicateId), path("${bam.baseName}.uniq.bam"), path("${bam.baseName}.uniq.bam.bai")
	
    script:
    """
    # Filter for unique alignments and create a new BAM file
    (samtools view -H $bam; samtools view $bam | grep -w 'NH:i:1') \
    | samtools view -Sb - > ${bam.baseName}.uniq.bam

    # Index the BAM file
    samtools index ${bam.baseName}.uniq.bam
    """
}

process ADD_READ_GROUP {
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/star/RG", mode: "copy"

    input:
    tuple val(replicateId), path(bamFile), path(baiFile) 

    output:
    tuple val(replicateId), path("rg_added.bam")

    script:
    """
    gatk AddOrReplaceReadGroups \
        -I $bamFile \
        -O rg_added.bam \
        -RGID '$replicateId' \
        -RGLB 'library' \
        -RGPL 'illumina' \
        -RGPU 'machine' \
        -RGSM 'GM12878'
    """
}



/*
 * Process 3: Split reads that contain Ns in their CIGAR string.
 *            Creates k+1 new reads (where k is the number of N cigar elements) 
 *            that correspond to the segments of the original read beside/between 
 *            the splicing events represented by the Ns in the original CIGAR.
 */

process RNASEQ_GATK_SPLITNCIGAR {
    tag "$replicateId"
    label 'mem_large'
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/SNG", mode: 'copy'

    input:
    path genome
    path index
    path genome_dict
    tuple val(replicateId), path(bam)

    output:
    tuple val(replicateId), path('split.bam')

    script:
    """
    gatk SplitNCigarReads \
        -R $genome \
        -I $bam \
        --refactor-cigar-string \
        -O split.bam
    """
}


// Process to run samtools indexing
process SAMTOOLS_INDEX_SPLIT_BAM {
    container "https://depot.galaxyproject.org/singularity/samtools%3A0.1.18--h50ea8bc_13"
    publishDir "${params.outdir}/SNG/splitbam", mode: 'copy'

    input:
    tuple val(replicateId), path(bam)

    output:
    tuple val(replicateId), path("split.bam"), path("split.bam.bai")

    script:
    """
    samtools index $bam
    """
}

/*
 * Process 4: Base recalibrate to detect systematic errors in base quality scores, 
 *            select unique alignments and index
 *             
 */

process RNASEQ_GATK_RECALIBRATE {
    tag "$replicateId"
    label "mem_large"
    container "/home/kothai/cq-git-sample/Praktikum/gatk4_4.2.6.0--hdfd78af_0.sif"
    publishDir "${params.outdir}/recalibrate", mode: 'copy'

    input:
    path genome
    path index
    path dict
    tuple val(replicateId), path(bam), path(bai)
    tuple path(variants_file), path(variants_file_index)

    output:
    tuple val(replicateId), path("${replicateId}.final.uniq.bam")

    script:
    """
    gatk BaseRecalibrator \
        -R $genome \
        -I $bam \
        --known-sites $variants_file \
        -O final.rnaseq.grp \
        --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true"

    gatk ApplyBQSR \
        -R $genome \
        -I $bam \
        --bqsr-recal-file final.rnaseq.grp \
        -O ${replicateId}.final.uniq.bam \
        --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
    """
}

process SAMTOOLS_INDEX_BAM {
    tag "$replicateId"
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.19.1--h50ea8bc_0"
    publishDir "${params.outdir}/recalibrate/bam", mode: 'copy'

    input:
    tuple val(replicateId), path(bam)

    output:
    tuple val(replicateId), path("${bam.baseName}.bam"), path("${bam.baseName}.bam.bai")

    script:
    """
    if [ -f $bam ]; then
        samtools index -b $bam ${bam.baseName}.bam.bai
    else
        echo "Error: BAM file not found" >&2
        exit 1
    fi
    """
}

/*
 * Process 5: Call variants with GATK HaplotypeCaller.
 *            Calls SNPs and indels simultaneously via local de-novo assembly of 
 *            haplotypes in an active region.
 *            Filter called variants with GATK VariantFiltration.    
 */

process RNASEQ_CALL_VARIANTS {
    tag "$sampleId"
    label "mem_xlarge"
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
        path genome         // Path to the genome reference (e.g., genome.fa)
        path genome_fai     // Path to the .fai index file for the genome
        path genome_dict    // Path to the genome dictionary file (e.g., genome.dict)
        tuple val(sampleId), path(bam), path(bai)

    output:
        tuple val(sampleId), path("final_${sampleId}.vcf")

    script:
    """
	
	# Print paths for debugging
    echo "Genome: ${genome}"
    echo "BAM file: ${bam}"
    echo "Output: output_${sampleId}.vcf.gz"
	
	
    # Adjust genome dictionary path to absolute if necessary
    sed -i 's|UR:file:.*|UR:file:/home/kothai/cq-git-sample/Praktikum/data/genome.fa|' ${genome_dict}

    # Variant calling with HaplotypeCaller
    gatk HaplotypeCaller \
        --native-pair-hmm-threads ${task.cpus} \
        --reference ${genome} \
        --output output_${sampleId}.vcf.gz \
        -I ${bam} \
        --standard-min-confidence-threshold-for-calling 20.0 \
        --dont-use-soft-clipped-bases

    # Variant filtering with VariantFiltration
    gatk VariantFiltration \
        -R ${genome} -V output_${sampleId}.vcf.gz \
        --cluster-window-size 35 --cluster-size 3 \
        --filter-name FS --filter-expression "FS > 30.0" \
        --filter-name QD --filter-expression "QD < 2.0" \
        -O final_${sampleId}.vcf
    """
}



process ANNOTATE_VARIANTS {
    tag "$sampleId"
    label "mem_large"
    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
        tuple val(sampleId), path(vcf)

    output:
        tuple val(sampleId), path("${sampleId}.annotated.vcf"), path("${sampleId}.summary.html")

    script:
    """
    java -Xmx16G -jar /mnt/snpEff/snpEff.jar \\
        -c ${params.snpeff_config} \\
        -v ${params.snpeff_db} ${vcf} > ${sampleId}.annotated.vcf

    java -Xmx16G -jar /mnt/snpEff/snpEff.jar \\
        -c ${params.snpeff_config} \\
        -v ${params.snpeff_db} ${vcf} -stats ${sampleId}.summary.html > /dev/null
    """
}



process multiqc {
    tag "MultiQC"
    publishDir "${params.outdir}/multiqc", mode: 'copy'
	container "https://depot.galaxyproject.org/singularity/multiqc%3A1.24.1--pyhdfd78af_0"
    input:
    path results_dir

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc ${results_dir} --outdir .
    """
}


workflow {
    // Step 1: Create a channel for paired-end reads
    reads_channel = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { sample_id, reads -> tuple(sample_id, reads.sort()) }

    // Step 2: Run FastQC on raw reads
    fastqc_raw_results = fastqc_raw(reads_channel)

    // Step 3: Run Fastp for trimming
    fastp_results = fastp(reads_channel)

    // Step 4: Extract and regroup trimmed reads
    fastp_trimmed_reads = fastp_results.trimmed_reads
        .map { sample_id, r1, r2 -> tuple(sample_id, [r1, r2]) }

    // Step 5: Run FastQC on trimmed reads
    fastqc_trimmed_results = fastqc_trimmed(fastp_trimmed_reads)

    // Step 6: Data Preparation - Prepare genome indices and filtered VCF
    def genome_index_samtools = PREPARE_GENOME_SAMTOOLS(params.genome)
    def genome_dict = PREPARE_GENOME_PICARD(params.genome)
    def star_genome_dir = PREPARE_STAR_GENOME_INDEX(params.genome)
    def filtered_vcf = PREPARE_VCF_FILE(params.variants, params.denylist)

    // Step 7: STAR RNA-Seq Mapping
    star_aligned_bam = RNASEQ_MAPPING_STAR(star_genome_dir, fastp_trimmed_reads)

    // Step 8: Generate alignment statistics
    alignment_stats = SAMTOOLS_FLAGSTAT(star_aligned_bam)

    // Step 9: Process mapped reads with Samtools, Add Read Group, and SplitNCigar
    filtered_bam = SAMTOOLS_FILTER_INDEX(star_aligned_bam)
    bam_with_rg = ADD_READ_GROUP(filtered_bam)
    split_bam = RNASEQ_GATK_SPLITNCIGAR(params.genome, genome_index_samtools, genome_dict, bam_with_rg)
    indexed_split_bam = SAMTOOLS_INDEX_SPLIT_BAM(split_bam)

    // Step 10: Base Recalibration and Indexing
    recalibrated_bam = RNASEQ_GATK_RECALIBRATE(params.genome, genome_index_samtools, genome_dict, indexed_split_bam, filtered_vcf)
    indexed_final_bam = SAMTOOLS_INDEX_BAM(recalibrated_bam)

    // Step 11: Variant Calling and Annotation
    variant_vcf = RNASEQ_CALL_VARIANTS(params.genome, genome_index_samtools, genome_dict, indexed_final_bam)
    annotated_vcf = ANNOTATE_VARIANTS(variant_vcf)

    // Step 12: Run MultiQC to aggregate all results
    multiqc_results = multiqc(Channel.fromPath("${params.outdir}"))
}






