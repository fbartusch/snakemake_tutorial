# build container
```
sudo singularity build snakemake_tutorial.simg Singularity
```

# shell into container
```
singularity shell --bind data/:/scif/data snakemake_tutorial.simg
```

# Map with bwa mem
```
mkdir /scif/data/mapped_reads
scif exec bwa bwa mem -o /scif/data/mapped_reads/r1_subset.sam /scif/data/index/ref.fa /scif/data/raw_reads/r1_subset.fq
```

# Sam -> Bam
```
scif exec samtools samtools view -Sb /scif/data/mapped_reads/r1_subset.sam > /scif/data/mapped_reads/r1_subset.bam
```
