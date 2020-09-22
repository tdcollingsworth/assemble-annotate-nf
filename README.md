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

## Example run command:
```
nextflow run unicycler_assemble-annotate.nf -c CONFIG/unicycler.config
```
