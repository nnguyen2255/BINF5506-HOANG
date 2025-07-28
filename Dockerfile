FROM continuumio/miniconda3

# Set up environment
RUN apt-get update && apt-get install -y curl git wget default-jre build-essential

# Install bioinformatics tools
RUN conda install -c bioconda -c conda-forge \
    snakemake bwa samtools gatk4 fastqc snpeff

# Set working directory
WORKDIR /pipeline

# Copy Snakemake files into image
COPY . /pipeline

# Set the command to run Snakemake pipeline
CMD ["snakemake", "--cores", "1"]
