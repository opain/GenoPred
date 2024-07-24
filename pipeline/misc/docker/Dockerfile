# Use the base image
FROM ubuntu:latest

# Run commands in a single RUN layer to reduce the image size
RUN apt-get update && \
    apt-get install -y git curl wget zip gzip vim build-essential zlib1g-dev jq && \

    # Create a wrapper script for tar that includes --no-same-owner
    echo '#!/bin/sh' > /usr/local/bin/tar && \
    echo 'exec /bin/tar --no-same-owner "$@"' >> /usr/local/bin/tar && \
    chmod +x /usr/local/bin/tar

# Install Mambaforge
RUN wget https://github.com/conda-forge/miniforge/releases/download/4.10.3-10/Mambaforge-4.10.3-10-Linux-x86_64.sh && \
    bash Mambaforge-4.10.3-10-Linux-x86_64.sh -b -p /opt/mambaforge && \
    rm Mambaforge-4.10.3-10-Linux-x86_64.sh

# Download GenoPred
RUN mkdir -p /tools/GenoPred && \
    cd /tools/GenoPred && \
    git clone https://github.com/opain/GenoPred.git . && \
    git fetch --tags

# Create genopred env
RUN /opt/mambaforge/bin/conda init bash && \
    echo "source /opt/mambaforge/etc/profile.d/conda.sh" >> /etc/bash.bashrc && \
    . /opt/mambaforge/etc/profile.d/conda.sh && \
    mamba env update -n genopred -f /tools/GenoPred/pipeline/envs/pipeline.yaml && \
    mamba clean -a -y

# Run initial snakemake command to setup dependencies
RUN /opt/mambaforge/bin/conda init bash && \
    echo "source /opt/mambaforge/etc/profile.d/conda.sh" >> /etc/bash.bashrc && \
    . /opt/mambaforge/etc/profile.d/conda.sh && \
    conda activate genopred && \
    cd /tools/GenoPred/pipeline && \
    snakemake --restart-times 3 -j1 --use-conda --conda-frontend mamba install_r_packages resources/software/pgscatalog_utils/download_pgscatalog_utils.done && \
    mamba clean -a -y

ENTRYPOINT ["/bin/bash"]
CMD ["-c", "source /opt/mambaforge/etc/profile.d/conda.sh && conda init && source ~/.bashrc && conda activate genopred && cd /tools/GenoPred/pipeline && exec /bin/bash"]
