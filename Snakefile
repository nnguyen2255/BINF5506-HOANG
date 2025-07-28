# Variables
SRA = "SRR1972739"
REF_ID = "AF086833.2"
RESULTS_FOLDER = "results"
RAW_DIR=f"{RESULTS_FOLDER}/raw"
ALIGNED_DIR=f"{RESULTS_FOLDER}/aligned"
VARIANT_DIR=f"{RESULTS_FOLDER}/variants"
ANNOTATED_DIR=f"{RESULTS_FOLDER}/annotated"
QC_DIR=f"{RESULTS_FOLDER}/qc"
SNPEFF_DIR=f"{RESULTS_FOLDER}/snpEff"
SNPEFF_DATA_DIR=f"{SNPEFF_DIR}/data/reference_db"
SNAKEMAKE_DIR=f"{RESULTS_FOLDER}/snakemake"

 
rule all:
    input: 
        f"{SNAKEMAKE_DIR}/.dirs_created",
        f"{RAW_DIR}/reference.fasta",
        f"{RAW_DIR}/{SRA}/{SRA}.sra",
        f"{RAW_DIR}/{SRA}.fastq",
        f"{QC_DIR}/{SRA}_fastqc.zip",
        f"{QC_DIR}/{SRA}_fastqc.html",
        f"{RAW_DIR}/reference.dict",
        f"{ALIGNED_DIR}/aligned.sam",
        f"{ALIGNED_DIR}/aligned.sorted.bam",
        f"{ALIGNED_DIR}/dedup.bam",
        f"{ALIGNED_DIR}/dup_metrics.txt",
        f"{VARIANT_DIR}/raw_variants.vcf",
        f"{VARIANT_DIR}/filtered_variants.vcf",
        f"{SNPEFF_DATA_DIR}/genes.gbk",
        f"{SNPEFF_DIR}/snpEff.config",
        f"{SNPEFF_DIR}/snpEff_reference_db.txt",
        f"{ANNOTATED_DIR}/annotated_variants.vcf"
rule create_dirs:
    output:
        marker = f"{SNAKEMAKE_DIR}/.dirs_created"
    shell:
        """
        mkdir -p {RESULTS_FOLDER} {RAW_DIR} {ALIGNED_DIR} {VARIANT_DIR} {ANNOTATED_DIR} {QC_DIR} {SNPEFF_DATA_DIR} {SNAKEMAKE_DIR}
        touch {output.marker}
        """
 
rule download_reference:
    input:
        marker = rules.create_dirs.output.marker
    output:
        reference_fasta = f"{RAW_DIR}/reference.fasta"
    shell:
        """
        echo Downloading FASTA reference...
        efetch -db nucleotide -id {REF_ID} -format fasta > {RAW_DIR}/reference.fasta
        echo Downloaded reference.fasta!
        """

 
rule download_sra:
    input:
        marker = rules.create_dirs.output.marker
    output:
        sequence_sra = f"{RAW_DIR}/{SRA}/{SRA}.sra"
    shell:
        """
        echo Downloading sequencing data...
        prefetch {SRA} -O {RAW_DIR}
        echo Downloaded sequencing data!
        """
 
rule extract_sequence:
    input:
        marker = rules.create_dirs.output.marker,
        sequence_sra = rules.download_sra.output.sequence_sra
    output:
        sequence_fastq = f"{RAW_DIR}/{SRA}.fastq"
    shell:
        """
        echo Extracting sequencing data...
        fastq-dump -X 10000 {RAW_DIR}/{SRA}/{SRA}.sra -O {RAW_DIR}
        echo Extracted sequencing data!
        """
 
rule run_fastQC:
    input:
        marker = rules.create_dirs.output.marker,
        sequence_sra = rules.download_sra.output.sequence_sra,
        extract_sequence = rules.extract_sequence.output.sequence_fastq
    output:
        zip = f"{QC_DIR}/{SRA}_fastqc.zip",
        html = f"{QC_DIR}/{SRA}_fastqc.html"
    shell:
        """
        echo Running FastQC on raw reads...
        fastqc -o {QC_DIR} {RAW_DIR}/{SRA}.fastq
        """

rule create_fasta_dict:
    input:
        marker = rules.create_dirs.output.marker,
        sequence_sra = rules.download_sra.output.sequence_sra,
        extract_sequence = rules.extract_sequence.output.sequence_fastq
    output:
        index = f"{RAW_DIR}/reference.dict"
    shell:
        f"""
        echo Indexing reference genome with samtools...
        samtools faidx {RAW_DIR}/reference.fasta

        echo Building BWA index...
        bwa index {RAW_DIR}/reference.fasta

        echo Creating FASTA dictionary using GATK...
        gatk CreateSequenceDictionary -R {RAW_DIR}/reference.fasta -O {RAW_DIR}/reference.dict
        """

rule alignment_read:
    input:
        marker = rules.create_dirs.output.marker,
        sequence_sra = rules.download_sra.output.sequence_sra,
        extract_sequence = rules.extract_sequence.output.sequence_fastq
    output:
        index = f"{ALIGNED_DIR}/aligned.sam"
    shell:
        f"""
        echo Aligning reads with read groups...
        bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:sample1' {RAW_DIR}/reference.fasta {RAW_DIR}/{SRA}.fastq > {ALIGNED_DIR}/aligned.sam
        echo Aligned reads!
        """
rule convert_sam_to_bam:
    input:
        index = rules.alignment_read.output.index
    output:
        bam_file = f"{ALIGNED_DIR}/aligned.sorted.bam"
    shell:
        """
        echo Converting SAM to sorted BAM...
        samtools view -b {ALIGNED_DIR}/aligned.sam | samtools sort -o {ALIGNED_DIR}/aligned.sorted.bam
        """

rule validate_duplicate_bam:
    input:
        bam_file = rules.convert_sam_to_bam.output.bam_file
    output:
        dedup_bam  = f"{ALIGNED_DIR}/dedup.bam",
        dup_metric = f"{ALIGNED_DIR}/dup_metrics.txt"
    shell:
        """
        echo Validating BAM file...
        gatk ValidateSamFile -I {ALIGNED_DIR}/aligned.sorted.bam -MODE SUMMARY
        
        echo Marking duplicates...
        gatk MarkDuplicates -I {ALIGNED_DIR}/aligned.sorted.bam -O {ALIGNED_DIR}/dedup.bam -M {ALIGNED_DIR}/dup_metrics.txt
        """

rule calling_variants:
    input:
        reference_fasta = rules.download_reference.output.reference_fasta,
        dedup_bam  = rules.validate_duplicate_bam.output.dedup_bam
    output:
        raw_variant = f"{VARIANT_DIR}/raw_variants.vcf"
    shell:
        """
        echo Indexing deduplicated BAM file...
        samtools index {ALIGNED_DIR}/dedup.bam

        echo Calling variants...
        gatk HaplotypeCaller -R {RAW_DIR}/reference.fasta -I {ALIGNED_DIR}/dedup.bam -O {VARIANT_DIR}/raw_variants.vcf
        echo Called variants!
        """

rule filter_variants:
    input:
        reference_fasta = rules.download_reference.output.reference_fasta,
        raw_variant = rules.calling_variants.output.raw_variant
    output:
        filter_variants = f"{VARIANT_DIR}/filtered_variants.vcf"

    shell:
        """
        echo Filtering variants...
        gatk VariantFiltration -R {RAW_DIR}/reference.fasta -V {VARIANT_DIR}/raw_variants.vcf -O {VARIANT_DIR}/filtered_variants.vcf --filter-expression "QD < 2.0 || FS > 60.0" --filter-name FILTER
        """

rule download_GenBank:
    input:
        marker = rules.create_dirs.output.marker
    output:
        genbank = f"{SNPEFF_DATA_DIR}/genes.gbk"

    shell:
        """
        echo Downloading reference GenBank file for snpEff...
        efetch -db nucleotide -id {REF_ID} -format genbank > {SNPEFF_DATA_DIR}/genes.gbk
        echo Downloaded GenBank file for snpEff!
        """

rule snpEff_config:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
        genbank = rules.download_GenBank.output.genbank
    output:
        snpEff_config = f"{SNPEFF_DIR}/snpEff.config"

    shell:
        """
        echo Creating custom snpEff configuration file...
        cat <<EOF > {SNPEFF_DIR}/snpEff.config
        # Custom snpEff config for reference_db
        reference_db.genome : reference_db
        reference_db.fa : $(readlink -f {RAW_DIR}/reference.fasta)
        reference_db.genbank : $(readlink -f {SNPEFF_DATA_DIR}/genes.gbk)
        EOF
        """

rule snpEff_database:
    input:
        snpEff_config = rules.snpEff_config.output.snpEff_config
    output:
        snpEff_database = f"{SNPEFF_DIR}/snpEff_reference_db.txt"

    shell:
        """
        echo Building snpEff database...
        snpEff build -c {SNPEFF_DIR}/snpEff.config -genbank -v -noCheckProtein reference_db

        echo Built snpEff database!
        echo Exporting snpEff database...
        snpEff dump -c {SNPEFF_DIR}/snpEff.config reference_db > {SNPEFF_DIR}/snpEff_reference_db.txt 
        echo Exported snpEff database!
        """

rule annotating_variants:
    input:
        reference_fasta = rules.download_reference.output.reference_fasta,
        filter_variants = rules.filter_variants.output.filter_variants,
        snpEff_config = rules.snpEff_config.output.snpEff_config,
        snpEff_database = rules.snpEff_database.output.snpEff_database,

    output:
        annotated_variants = f"{ANNOTATED_DIR}/annotated_variants.vcf"

    shell:
        """
        echo Annotating variants with snpEff...
        snpEff -c {SNPEFF_DIR}/snpEff.config -stats {SNPEFF_DIR}/snpEff.html reference_db {VARIANT_DIR}/filtered_variants.vcf > {ANNOTATED_DIR}/annotated_variants.vcf
        echo Annotated variants with snpEff!

        echo Pipeline completed successfully! Check the folders in {RESULTS_FOLDER} for output files.
        tree
        echo All done!
        """