#variables:
SRA=SRR1972739
REF_ID=AF086833.2
RESULTS_FOLDER=results
RAW_DIR=${RESULTS_FOLDER}/raw
ALIGNED_DIR=${RESULTS_FOLDER}/aligned
VARIANT_DIR=${RESULTS_FOLDER}/variants
ANNOTATED_DIR=${RESULTS_FOLDER}/annotated
QC_DIR=${RESULTS_FOLDER}/qc
SNPEFF_DIR=${RESULTS_FOLDER}/snpEff
SNPEFF_DATA_DIR=${SNPEFF_DIR}/data/reference_db

rule all:
    input:
        "hello.txt",
        "goodbye.txt"

rule hello: 
    output:
        "hello.txt"
    shell:
        "echo Hello World" > {output}

rule goodbye:
    output:
        "goodbye.txt"
    shell:
        "echo Goodbye World" > {output}