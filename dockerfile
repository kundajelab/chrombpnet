# Use the official TensorFlow image as parent
	FROM tensorflow/tensorflow:2.5.1-gpu

# Set the working directory
WORKDIR /scratch


# Install some basic utilities
RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git && \
	    apt-get clean && \
		    rm -rf /var/lib/apt/lists/*

# Install Google Cloud SDK
RUN apt-get update && apt install -y --allow-unauthenticated wget
RUN cd /opt/ && \
	wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-307.0.0-linux-x86_64.tar.gz && \
		tar xvfz google-cloud-sdk-307.0.0-linux-x86_64.tar.gz && \
			./google-cloud-sdk/install.sh
			ENV PATH "$PATH:/opt/google-cloud-sdk/bin/"

# Install Miniconda with Python 3.7 into /opt
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_4.9.2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
	    rm ~/miniconda.sh && \
		    /opt/conda/bin/conda clean -tipsy

# Enable Conda and alter bashrc so the Conda default environment is always activated
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
	    echo "conda activate base" >> ~/.bashrc

# Attach Conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Copy the entire repo
RUN mkdir /scratch/chrombpnet-lite
COPY . /scratch/chrombpnet-lite

# Install any needed packages specified in requirements.txt
RUN pip install -r /scratch/chrombpnet-lite/requirements.txt

# Set environment variables for Python
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

