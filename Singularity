Bootstrap: docker
From: continuumio/anaconda3

# sudo singularity build hello-world-scif.simg Singularity.scif

%files
    snakemake_tutorial.scif

%environment
    PATH=/opt/conda/bin:$PATH
    export PATH

%post
    apt-get -y install build-essential

    # Install scif, snakemake
    /opt/conda/bin/pip install scif 
    /opt/conda/bin/pip install snakemake==4.4.0
    /opt/conda/bin/scif install snakemake_tutorial.scif

    mkdir /input
    mkdir /output

%runscript
    exec scif "$@"
