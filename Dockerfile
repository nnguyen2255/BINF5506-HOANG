FROM continuumio/miniconda3

# Set up working directory
WORKDIR /pipeline

# Copy the conda environment file (see next step)
COPY environment.yml /tmp/environment.yml

# Create the conda environment
RUN conda env create -f /tmp/environment.yml

# Install EDirect (efetch)
RUN apt-get update && \
    apt-get install -y perl wget unzip && \
    mkdir -p /root/edirect && \
    cd /root/edirect && \
    wget https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz && \
    tar -xzvf edirect.tar.gz && \
    rm edirect.tar.gz && \
    echo 'export PATH=$PATH:/root/edirect' >> /root/.bashrc

# Ensure PATH includes edirect tools at runtime
ENV PATH="/root/edirect:$PATH"

# Copy your pipeline files into the image
COPY . /pipeline

# Run Snakemake as the default command
CMD ["conda", "run", "--no-capture-output", "-n", "snakemake_env", "snakemake", "--cores", "1"]
