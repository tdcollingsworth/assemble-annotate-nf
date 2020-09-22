#!/usr/bin/env nextflow

// Script parameters
params.forward_reads = "${launchDir}/DATA/atcc2001_S12_L001_R1_001.fastq.gz"
params.reverse_reads = "${launchDir}/DATA/atcc2001_S12_L001_R2_001.fastq.gz"
params.long_reads = "${launchDir}/DATA/nanopore.fasta"
params.reference = "${launchDir}/DATA/GCF_000002545.3_ASM254v2_protein.faa" 

forward_ch = Channel.fromPath(params.forward_reads)
reverse_ch = Channel.fromPath(params.reverse_reads)
nanopore_ch = Channel.fromPath(params.long_reads)
reference_ch = Channel.fromPath(params.reference)

// Processes

process masurca {

    label 'porky'

    input:
    file forward from forward_ch
    file reverse from reverse_ch
    file nanopore from nanopore_ch

    output:
    file "assembly.fasta" into assembly_ch

    """
    # Output directroy
    mkdir MASURCA ; cd MASURCA/
    # Get configuration script, edit
    cp /TOOLS/masurca/sr_config_example.txt config.txt
    sed -i "s|PE= pe 500 50  /FULL_PATH/frag_1.fastq  /FULL_PATH/frag_2.fastq|PE= pe 500 50 ${launchDir}/DATA/${forward} ${launchDir}/DATA/${reverse}|g" config.txt
    sed -i "s|JUMP= sh 3600 200 |#JUMP= sh 3600 200 |g" config.txt
    sed -i "s|#NANOPORE=/FULL_PATH/nanopore.fa|NANOPORE=${launchDir}/DATA/${nanopore}|g" config.txt
    sed -i "s|EXTEND_JUMP_READS=0|#EXTEND_JUMP_READS=0|g" config.txt
    sed -i "s|FLYE_ASSEMBLY=0|FLYE_ASSEMBLY=1|g" config.txt
    cat config.txt
    # Set up masurca and run    
    masurca config.txt
    ./assemble.sh
    """
}

process annotate {

    label 'porky'

    input:
    file assembly from assembly_ch
    file reference from reference_ch

    output:
    file "augustus.hints.gtf" into annotation_ch

    """
    ## Simplify reference '.faa' headers ...
    simplifyFastaHeaders.pl ${reference} protein_aa_ protein_aa.fa protein_aa_header.map
    
    ## Soft mask the draft genome ...
    BuildDatabase -engine ncbi -name assembly ${assembly}
    RepeatModeler -database assembly -engine ncbi -pa ${task.cpus}
    RepeatMasker -pa ${task.cpus} -e ncbi -s -xsmall -lib assembly-families.fa ${assembly}
    
    ## Simplify masked assembly '.fasta' headers ...
    simplifyFastaHeaders.pl assembly.fasta.masked assembly_ assembly.fa assembly_header.map
    
    ## Run BRAKER
    braker.pl --genome=/WORKSPACE/assembly.fa --prot_seq=/WORKSPACE/protein_aa.fa --prg=gth --trainFromGth --softmasking --species=assembly --verbosity=4 --workingdir=/WORKSPACE/ANNOTATE/ --nocleanup --cores=${task.cpus}
    """
}

