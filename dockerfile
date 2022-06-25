# Use the official TensorFlow image as parent
FROM tensorflow/tensorflow:2.8.2-gpu

# Set the working directory
WORKDIR /scratch

# Install some basic utilities
RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git

# Install Google Cloud SDK
RUN apt-get update && apt install -y --allow-unauthenticated wget
RUN cd /opt/ && \
    wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-307.0.0-linux-x86_64.tar.gz && \
    tar xvfz google-cloud-sdk-307.0.0-linux-x86_64.tar.gz && \
    ./google-cloud-sdk/install.sh
ENV PATH "$PATH:/opt/google-cloud-sdk/bin/"

# Install Miniconda with Python 3.7 into /opt
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Enable Conda and alter bashrc so the Conda default environment is always activated
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Attach Conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Install SAMtools, BEDtools, and UCSC BedGraphToBigWig
RUN conda install -y -c bioconda samtools bedtools ucsc-bedgraphtobigwig 

# Install jq
RUN apt-get install -y jq

# Clean up after apt and conda
RUN apt-get clean && rm -rf /var/lib/apt/lists/*
RUN conda clean -tipy

# Set environment variables for Python
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# Copy the entire repo
RUN mkdir /scratch/chrombpnet
COPY . /scratch/chrombpnet

# Install any needed packages specified in requirements.txt
RUN pip install -r /scratch/chrombpnet/requirements.txt


