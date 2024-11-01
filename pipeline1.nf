nextflow.enable.dsl = 2
//ALL INPUT DATA 
//1.REFERENCE GENOME ASSEMBLY(.FA)
//2.RNA SEQ READS(FASTAQ.GZ)
//3. KNOWN VARIANTS(VCF.GZ)
//4.BLACKLISTED REGIONS (.BED)
//LOT OF INTERMEDIATE FILES AND OUTPUTS

//MAIN AIM OF THIS PIPELINE IS TO PROCESS THE RAW RNA-SEQ DATA AND OBTAIN THE LIST OF SMALL VARIANTS FOR DOWNSTREAM ANALYSIS. 
//PIPELINE INCLUDES SNVS POSTPROCESSING AND QUANTIFICATION FOR ALLELE SPECIFIC EXPRESSION. 

//FIRST STEP IS TO PREPARE THE GENOME
//GENOME INDICES WITH SAMTOOLS AND PICARD IS DONE FIRST. THEW WILL BE NEEDED FOR GATK COMMANDS SUCH AS SPLITNTRIM
//GENOME INDEX FOR STAR IS CREATED NEXTFLOW
//VARIANT OVERLAPPING BLACKLISTED REGIONS ARE FILTERED THEN TO REDUCE FALSE POSITIVE CALLS

//MAPPING (2-PASS APPROACH). FIRST ALIGNMENT CREATES TABLE WITH SPLICE JUNCTIONS THAT IS USED TO GUIDE FINAL ALIGNMENT. 
//BAM FILES - CREATE INDEX. 

//BASE RECALLIBRATION
//VARIANT CALLING AND VARIANT FILTERING
//VARIANT POST PROCESSING. 


//STEP1: 
// DEFINING ALL THE PARAMETERS
params.genome = "$launchDir/data/genome.fa"
params.variants = "$launchDir/data/known_variants.vcf.gz"
params.denylist   = "$launchDir/data/denylist.bed"
params.results = "results"
//params.transcriptome_file = "$launchDir/data/ggal/transcriptome.fa"
params.accession = ""
params.with_fastqc = true
params.with_fastp = true
params.cache = "${launchDir}/cache"
params.out = "${launchDir}/out"
params.reads = "$baseDir/data/reads/rep1_{1,2}.fq.gz"

process prefetch {
  container "https://depot.galaxyproject.org/singularity/sra-tools%3A2.11.0--pl5321ha49a11a_3"
  storeDir params.cache
  input:
    val accession
  output:
    path "${accession}"
  """
  prefetch ${accession}
  """
}

process fasterqdump {
  container "https://depot.galaxyproject.org/singularity/sra-tools%3A2.11.0--pl5321ha49a11a_3"
  storeDir params.cache
  input:
    path sradir
  output:
    path "${sradir}_*.fastq"
  """
  fasterq-dump --split-files ${sradir} 
  """
} 


// Step 3: Run Fastp (Trimming)
process fastp {
	storeDir params.cache
	publishDir params.out, mode: 'copy', overwrite: true
	container "https://depot.galaxyproject.org/singularity/fastp%3A0.23.4--hadf994f_3"
	input:
		tuple val(sampleId), path(reads) 
	output:
		tuple val(sampleId), path("trimmed_*.fastq.gz")
	
	script:
	"""
	echo "Processing sample ID: ${sampleId}"
    echo "Input FASTQ files: ${reads[0]} and ${reads[1]}"
	
	mkdir -p ${params.out}/fastp
    fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 trimmed_1.fastq.gz --out2 trimmed_2.fastq.gz \
          --html ${params.out}/fastp/fastp_report_${sampleId}.html --json ${params.out}/fastp/fastp_report_${sampleId}.json
	"""
}

// Step 5: Run FastQC on trimmed reads
process fastqc_trimmed {
	publishDir params.out, mode: "copy", overwrite: true
	container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
	input:
        tuple val(sampleId), path(trimmed_fastq_files)
    
    output:
        tuple val(sampleId), path("*_fastqc.zip"), path("*_fastqc.html") 
    
    script:
    """
	echo "Running FastQC on sample ID: ${sampleId}"
    echo "Trimmed FASTQ files: ${trimmed_fastq_files[0]} and ${trimmed_fastq_files[1]}"

	fastqc ${trimmed_fastq_files[0]} ${trimmed_fastq_files[1]} 
	"""
}
//Create a fasta Genome Index

process PREPARE_GENOME_SAMTOOLS {
	container "https://depot.galaxyproject.org/singularity/samtools%3A0.1.18--h50ea8bc_13"
	publishDir "data", mode: 'copy', overwrite: true
	
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
	publishDir "data", mode: 'copy', overwrite: true
	
	input:
		path genome
	output:
		path "${genome.baseName}.dict"

	script:
	"""
	gatk CreateSequenceDictionary \
		-R $genome \
		-O ${genome.baseName}.dict 
    
	# Correct any UR:file path in the resulting genome.dict file to absolute path
    sed -i 's|UR:file:.*|UR:file:/home/kothai/cq-git-sample/Praktikum/data/genome.fa|' ${genome.baseName}.dict
	"""
}

//Create a genome Index file

process PREPARE_STAR_GENOME_INDEX {
	container "https://depot.galaxyproject.org/singularity/star%3A2.7.10b--h6b7c446_1"
	publishDir params.out, mode: 'copy', overwrite: true
	
	tag "$genome.baseName"
	input:
		path genome
	output:
		path "genome_dir"
	script:
	"""
	mkdir genome_dir 
	
	STAR --runMode genomeGenerate \
		 --genomeDir genome_dir \
		 --genomeFastaFiles ${genome} \
		 --runThreadN ${task.cpus}
	"""
}

//Filtered and recorded variants

process PREPARE_VCF_FILE {
    container "https://depot.galaxyproject.org/singularity/mulled-v2-b9358559e3ae3b9d7d8dbf1f401ae1fcaf757de3%3Aac05763cf181a5070c2fdb9bb5461f8d08f7b93b-0"
     publishDir params.out, mode: 'copy', overwrite: true
	 
    input: 
		path variantsFile
		path denylisted

	output:
    tuple \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz"), \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")
  
	script:  
	"""
	vcftools --gzvcf $variantsFile -c \
           --exclude-bed ${denylisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz

	tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
	"""
}

process RNASEQ_MAPPING_STAR {
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
    
    publishDir params.out, mode: 'copy', overwrite: true
	
	tag "$replicateId"

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


process SAMTOOLS_FILTER_INDEX {
    container "https://depot.galaxyproject.org/singularity/samtools%3A0.1.18--h50ea8bc_13"
    

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
    container "/home/kothai/cq-git-sample/Praktikum/gatk4_4.2.6.0--hdfd78af_0.sif"

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
	
    publishDir params.out, mode: 'copy', overwrite: true
	container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
	input: 
		path genome
		path index
		path genome_dict
		tuple val(replicateId), path(bam)


	output:
		tuple val(replicateId), path('split.bam')
  
	script:
	"""
	# SplitNCigarReads and reassign mapping qualities
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
    publishDir params.out, mode: 'copy', overwrite: true

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
    publishDir params.out, mode: 'copy', overwrite: true
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
    # Indel Realignment and Base Recalibration
    gatk BaseRecalibrator \
          -R $genome \
          -I $bam \
          --known-sites $variants_file \
          -O final.rnaseq.grp \
          --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true"

    gatk ApplyBQSR \
          -R $genome -I $bam \
          --bqsr-recal-file final.rnaseq.grp \
          -O ${replicateId}.final.uniq.bam \
          --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
    """
}
process SAMTOOLS_INDEX_BAM {
    tag "$replicateId"
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.19.1--h50ea8bc_0"
	publishDir params.out, mode: 'copy', overwrite: true
	
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
    container "/home/kothai/cq-git-sample/Praktikum/gatk4_4.2.6.0--hdfd78af_0.sif"
    publishDir params.out, mode: 'copy', overwrite: true

    input:
        path genome         // Path to the genome reference (e.g., genome.fa)
        path genome_fai     // Path to the .fai index file for the genome
        path genome_dict    // Path to the genome dictionary file (e.g., genome.dict)
        tuple val(sampleId), path(bam)

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
    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.0--hdfd78af_1"
	publishDir params.out, mode: 'copy', overwrite: true
    input:
        tuple val(sampleId), path(vcf)

    output:
        tuple val(sampleId), path("${sampleId}.annotated.vcf")

    script:
    """
    snpEff -Xmx16G -c /home/kothai/cq-git-sample/Praktikum/snpEff/snpEff.config -v GRCh38.86 ${vcf} > ${sampleId}.annotated.vcf
    """
}




	
// Workflow definition
workflow {
    // Step 1: Define reads_channel based on whether params.accession is provided
    def reads_channel
    if (params.accession) {
        // Use prefetch and fasterqdump if accession is provided
        def accession_channel = Channel.value(params.accession)
        def sra_channel = prefetch(accession_channel)
        reads_channel = fasterqdump(sra_channel)
    } else {
        // Directly use params.reads if no accession is provided
        reads_channel = Channel.fromFilePairs(params.reads, checkIfExists: true)
    }

    // Step 2: Pass reads_channel to fastp
    def trimmed_reads = fastp(reads_channel)
    
    // Step 3: Optionally pass trimmed reads to FastQC if enabled
    if (params.with_fastqc) {
        fastqc_trimmed(trimmed_reads)
    }

    // Part 1: Data Preparation - Prepare genome indices and filtered VCF
    def genome_index_samtools = PREPARE_GENOME_SAMTOOLS(params.genome)
    def genome_dict = PREPARE_GENOME_PICARD(params.genome)
    def star_genome_dir = PREPARE_STAR_GENOME_INDEX(params.genome)
    def filtered_vcf = PREPARE_VCF_FILE(params.variants, params.denylist)

    // Part 2: STAR RNA-Seq Mapping
    def aligned_bam = RNASEQ_MAPPING_STAR(star_genome_dir, trimmed_reads)

    // Part 3: Process mapped reads with Samtools, add read group, and SplitNCigar
    def filtered_bam = SAMTOOLS_FILTER_INDEX(aligned_bam)
    def bam_with_rg = ADD_READ_GROUP(filtered_bam)
    def split_bam_output = RNASEQ_GATK_SPLITNCIGAR(params.genome, genome_index_samtools, genome_dict, bam_with_rg)
    def indexed_split_bam = SAMTOOLS_INDEX_SPLIT_BAM(split_bam_output)

    // Part 4: Base Recalibration and Indexing
    def recalibrated_bam = RNASEQ_GATK_RECALIBRATE(params.genome, genome_index_samtools, genome_dict, indexed_split_bam, filtered_vcf)
    def indexed_bam_final = SAMTOOLS_INDEX_BAM(recalibrated_bam)

    // Part 5: Variant Calling and Annotation
    def vcf_output = RNASEQ_CALL_VARIANTS(params.genome, genome_index_samtools, genome_dict, indexed_bam_final)
    def annotated_vcf = ANNOTATE_VARIANTS(vcf_output)
}






