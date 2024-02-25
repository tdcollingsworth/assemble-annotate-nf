:warning: In Development

# Unicycler
https://github.com/rrwick/Unicycler
```
docker build --no-cache -t collingswortht/unicycler -f DOCKER/UNICYCLER/Dockerfile-unicycler .
```

# Smartdenovo
https://github.com/ruanjue/smartdenovo
```
docker build --no-cache -t collingswortht/smartdenovo -f DOCKER/SMARTDENOVO/Dockerfile-smartdenovo .
```

# MaSuRCA
https://github.com/alekseyzimin/masurca/
```
docker build --no-cache -t collingswortht/masurca -f DOCKER/MASURCA/Dockerfile-masurca .
```

# BRAKER
https://github.com/Gaius-Augustus/BRAKER
```
docker build --no-cache -t collingswortht/annotate -f DOCKER/ANNOTATE/annotate_dockerfile .
```

## Example run:

Collect read data (short reads, forward and reverse, as well as long reads "\*.fastq.gz") and a reference fasta ("\*.faa").

Place in "DATA/" directory and edit "\*.nf" parameters.

```
nextflow run masurca_assemble-annotate.nf -c CONFIG/masurca.config
```
