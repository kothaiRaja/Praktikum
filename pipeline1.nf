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
params.accession = "SRR062634"
params.with_fastqc = false
params.with_fastp = false
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
		path fastq_files
	output:
		path "trimmed_*.fastq.gz"
	
	script:
	"""
	mkdir -p ${params.out}/fastp
	fastp --in1 ${fastq_files[0]} --in2 ${fastq_files[1]} --out1 trimmed_1.fastq.gz --out2 trimmed_2.fastq.gz \\
			--html ${params.out}/fastp/fastp_report.html --json ${params.out}/fastp/fastp_report.json
	"""
}

// Step 5: Run FastQC on trimmed reads
process fastqc_trimmed {
	publishDir params.out, mode: "copy", overwrite: true
	container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
	input:
        path trimmed_fastq_files
    
    output:
        path "*"  // FastQC output for reports
    
    script:
    """

	fastqc ${trimmed_fastq_files[0]} ${trimmed_fastq_files[1]} 
	"""
}
//Create a fasta Genome Index

process PREPARE_GENOME_SAMTOOLS {
	container "https://depot.galaxyproject.org/singularity/samtools%3A0.1.18--h50ea8bc_13"
	publishDir params.out, mode: 'copy', overwrite: true
	storeDir params.cache
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
	publishDir params.out, mode: 'copy', overwrite: true
	storeDir params.cache
	input:
		path genome
	output:
		path "${genome.baseName}.dict"

	script:
	"""
	gatk CreateSequenceDictionary -R $genome -O ${genome.baseName}.dict
	"""
}

//Create a genome Index file

process PREPARE_STAR_GENOME_INDEX {
	container "https://depot.galaxyproject.org/singularity/star%3A2.7.10b--h6b7c446_1"
	publishDir params.out, mode: 'copy', overwrite: true
	storeDir params.cache
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
	 storeDir params.cache
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
        tuple val(replicateId), path("${bam.baseName}.bai")

    script:
    """
    if [ -f $bam ]; then
        samtools index -b $bam ${bam.baseName}.bai
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

	input:
		path genome
		path index
		path dict
		tuple val(sampleId), path(bam), path(bai)
 
	output: 
		tuple val(sampleId), path('final.vcf')

	script:
	def bam_params = bam.collect{ "-I $it" }.join(' ')
	"""
	# fix absolute path in dict file
	sed -i 's@UR:file:.*${genome}@UR:file:${genome}@g' $dict
  
	# Variant calling
	gatk HaplotypeCaller \
          --native-pair-hmm-threads ${task.cpus} \
          --reference ${genome} \
          --output output.gatk.vcf.gz \
          ${bam_params} \
          --standard-min-confidence-threshold-for-calling 20.0 \
          --dont-use-soft-clipped-bases 

	# Variant filtering
	gatk VariantFiltration \
          -R ${genome} -V output.gatk.vcf.gz \
          --cluster-window-size 35 --cluster-size 3 \
          --filter-name FS --filter-expression \"FS > 30.0\" \
          --filter-name QD --filter-expression \"QD < 2.0\" \
          -O final.vcf
	"""
}

// Variant Annotation with SnpEff
process annotate_variants {
    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.0--hdfd78af_0"
    publishDir params.out, mode: 'copy', overwrite: true

    input:
        path vcf_file

    output:
        path "annotated_variants.vcf"

    script:
    """
    # Replace 'GRCh38.99' with the reference genome version you're using
    snpEff GRCh38.99 ${vcf_file} > annotated_variants.vcf
    """
}


/*
 * Process 6A: Post-process the VCF result  
 */

process POST_PROCESS_VCF {
	tag "$sampleId"
	publishDir "$params.results/$sampleId" 

	input:
		tuple val(sampleId), path('final.vcf')
		tuple path('filtered.recode.vcf.gz'), path('filtered.recode.vcf.gz.tbi')
	output: 
		tuple val(sampleId), path('final.vcf'), path('commonSNPs.diff.sites_in_files')
  
	script:
	'''
	grep -v '#' final.vcf | awk '$7~/PASS/' |perl -ne 'chomp($_); ($dp)=$_=~/DP\\=(\\d+)\\;/; if($dp>=8){print $_."\\n"};' > result.DP8.vcf
  
	vcftools --vcf result.DP8.vcf --gzdiff filtered.recode.vcf.gz  --diff-site --out commonSNPs
	'''
}

/* 
 * Process 6B: Prepare variants file for allele specific expression (ASE) analysis
 */

process PREPARE_VCF_FOR_ASE {
	tag "$sampleId"
	publishDir "$params.results/$sampleId" 

	input: 
		tuple val(sampleId), path('final.vcf'), path('commonSNPs.diff.sites_in_files')
	output: 
		tuple val(sampleId), path('known_snps.vcf.gz'), path('known_snps.vcf.gz.tbi')
		path 'AF.histogram.pdf'

	script:
	'''
	awk 'BEGIN{OFS="\t"} $4~/B/{print $1,$2,$3}' commonSNPs.diff.sites_in_files  > test.bed
    
	vcftools --vcf final.vcf --bed test.bed --recode --keep-INFO-all --stdout > known_snps.vcf

	grep -v '#'  known_snps.vcf | awk -F '\\t' '{print $10}' \
               |awk -F ':' '{print $2}'|perl -ne 'chomp($_); \
               @v=split(/\\,/,$_); if($v[0]!=0 ||$v[1] !=0)\
               {print  $v[1]/($v[1]+$v[0])."\\n"; }' |awk '$1!=1' \
               >AF.4R

	gghist.R -i AF.4R -o AF.histogram.pdf
	# Known SNPs have to be zipped and indexed for being used
	bgzip -c known_snps.vcf  > known_snps.vcf.gz
	tabix -p vcf known_snps.vcf.gz
	'''
}

/* 
 * Process 6C: Allele-Specific Expression analysis with GATK ASEReadCounter.
 *             Calculates allele counts at a set of positions after applying 
 *             filters that are tuned for enabling allele-specific expression 
 *             (ASE) analysis
 */

process ASE_KNOWNSNPS {
	tag "$sampleId"
	publishDir "$params.results/$sampleId" 
	label "mem_large"

	input:
		path genome
		path index
		path dict
		tuple val(sampleId), path(vcf), path(tbi), path(bam), path(bai)
  
	output:
		path "ASE.tsv"
  
	script:
	def bam_params = bam.collect{ "-I $it" }.join(' ')
	"""
	gatk ASEReadCounter \
          -R ${genome} \
          -O ASE.tsv \
          ${bam_params} \
          -V ${vcf}
	"""
}

	
// Workflow definition
workflow {
    // Step 1: Prefetch SRA file and convert to FASTQ format
    accession_channel = Channel.from(params.accession)
    sra_channel = prefetch(accession_channel)
    reads_channel = fasterqdump(sra_channel)

    // Optional Step 2: Trimming using Fastp (only if params.with_fastp is true)
    fastp_result = reads_channel
    if (params.with_fastp) {
        fastp_result = fastp(reads_channel)
    }

    // Optional Step 3: Run FastQC on trimmed reads (only if params.with_fastqc is true)
    if (params.with_fastqc) {
        fastqc_trimmed(fastp_result)
    }

    reads_ch = Channel.fromFilePairs(params.reads)

      // PART 1: Data preparation
      genome_index_samtools = PREPARE_GENOME_SAMTOOLS(params.genome)
      genome_dict = PREPARE_GENOME_PICARD(params.genome)
      star_genome_dir = PREPARE_STAR_GENOME_INDEX(params.genome)
      filtered_vcf = PREPARE_VCF_FILE(params.variants, params.denylist)
	  
	  
	  replicate_id_ch = Channel.value("rep1")

      // PART 2: STAR RNA-Seq Mapping
     aligned_bam = RNASEQ_MAPPING_STAR(  
            star_genome_dir, 
            reads_ch )
	
	 // Run Samtools filtering and indexing
     indexed_BAM = SAMTOOLS_FILTER_INDEX(aligned_bam)
	 
	 
	 
	 // Add read group information
    bam_with_rg = ADD_READ_GROUP(indexed_BAM)

      // PART 3: GATK Prepare Mapped Reads
      split_bam_output = RNASEQ_GATK_SPLITNCIGAR(
            params.genome, 
            PREPARE_GENOME_SAMTOOLS.out, 
            PREPARE_GENOME_PICARD.out, 
            bam_with_rg)
	 indexed_split_bam = SAMTOOLS_INDEX_SPLIT_BAM(split_bam_output)
	
	
      // PART 4: GATK Base Quality Score Recalibration Workflow
      recalibrated_bam = RNASEQ_GATK_RECALIBRATE(
                  params.genome, 
                  PREPARE_GENOME_SAMTOOLS.out, 
                  PREPARE_GENOME_PICARD.out, 
                  indexed_split_bam, 
                  PREPARE_VCF_FILE.out)
	 // Indexing the recalibrated BAM file
    indexed_bam = SAMTOOLS_INDEX_BAM(recalibrated_bam)

      // PART 5: GATK Variant Calling
      //RNASEQ_CALL_VARIANTS( 
            //params.genome, 
            //PREPARE_GENOME_SAMTOOLS.out, 
            //PREPARE_GENOME_PICARD.out, 
           // RNASEQ_GATK_RECALIBRATE.out.groupTuple() )

      // PART 6: Post-process variants file and prepare for 
      // Allele-Specific Expression and RNA Editing Analysis
      //POST_PROCESS_VCF( 
            //RNASEQ_CALL_VARIANTS.out, 
            //PREPARE_VCF_FILE.out )

      //PREPARE_VCF_FOR_ASE( POST_PROCESS_VCF.out )

      //ASE_KNOWNSNPS(
            //params.genome, 
            //PREPARE_GENOME_SAMTOOLS.out, 
            //PREPARE_GENOME_PICARD.out, 
            //group_per_sample(
                  //RNASEQ_GATK_RECALIBRATE.out, 
                  //PREPARE_VCF_FOR_ASE.out[0]) )
}




