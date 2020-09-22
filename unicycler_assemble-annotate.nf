#!/usr/bin/env nextflow

// Script parameters
params.forward_reads = "${launchDir}/DATA/atcc2001_S12_L001_R1_001.fastq.gz"
params.reverse_reads = "${launchDir}/DATA/atcc2001_S12_L001_R2_001.fastq.gz"
params.long_reads = "${launchDir}/DATA/atcc2001_nanopore.fastq.gz"
params.reference = "${launchDir}/DATA/GCF_000002545.3_ASM254v2_protein.faa" 

forward_ch = Channel.fromPath(params.forward_reads)
reverse_ch = Channel.fromPath(params.reverse_reads)
nanopore_ch = Channel.fromPath(params.long_reads)
reference_ch = Channel.fromPath(params.reference)

// Processes

process unicycler {

    label 'porky'

    input:
    file forward from forward_ch
    file reverse from reverse_ch
    file nanopore from nanopore_ch

    output:
    file "assembly.fasta" into assembly_ch

    """
    unicycler -t ${task.cpus} -o ${launchDir}/UNICYCLYER -1 ${forward} -2 ${reverse} -l ${nanopore}
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
    # Output directroy
    mkdir BRAKER ; cd BRAKER/

    # Simplify reference '.faa' headers ...
    simplifyFastaHeaders.pl ${reference} protein_aa_ protein_aa.fa protein_aa_header.map
    
    # Soft mask the draft genome ...
    BuildDatabase -engine ncbi -name assembly ${assembly}
    RepeatModeler -database assembly -engine ncbi -pa ${task.cpus}
    RepeatMasker -pa ${task.cpus} -e ncbi -s -xsmall -lib assembly-families.fa ${assembly}
    
    # Simplify masked assembly '.fasta' headers ...
    simplifyFastaHeaders.pl assembly.fasta.masked assembly_ assembly.fa assembly_header.map
    
    # Run BRAKER
    braker.pl --genome=/WORKSPACE/assembly.fa --prot_seq=/WORKSPACE/protein_aa.fa --prg=gth --trainFromGth --softmasking --species=assembly --verbosity=4 --workingdir=/WORKSPACE/ANNOTATE/ --nocleanup --cores=${task.cpus}
    """
}

