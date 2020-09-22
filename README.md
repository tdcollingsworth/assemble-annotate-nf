:warning: Still in early development

# Unicycler
https://github.com/rrwick/Unicycler
```
docker build --no-cache -t collingswortht/dev_unicycler -f DOCKER/UNICYCLER/Dockerfile-unicycler .
```

# Smartdenovo
https://github.com/ruanjue/smartdenovo
```
docker build --no-cache -t collingswortht/dev_smartdenovo -f DOCKER/SMARTDENOVO/Dockerfile-smartdenovo .
```

# MaSuRCA
https://github.com/alekseyzimin/masurca/
```
docker build --no-cache -t collingswortht/dev_masurca -f DOCKER/MASURCA/Dockerfile-masurca .
```

# BRAKER
https://github.com/Gaius-Augustus/BRAKER
```
docker build --no-cache -t collingswortht/dev_annotate -f DOCKER/ANNOTATE/annotate_dockerfile .
```

## Example run:

Collect read data (forward and reverse short reads as well as long reads "\*.fastq.gz") and a reference fasta ("\*.faa").

Place in "DATA/" directory and edit "\*.nf" parameters.

```
nextflow run unicycler_assemble-annotate.nf -c CONFIG/unicycler.config --
```
