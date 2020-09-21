#!/usr/bin/env nextflow

// Script parameters

params.forward_reads = "${launchDir}/DATA/atcc2001_S12_L001_R1_001.fastq.gz"
forward_ch = Channel.fromPath(params.forward_reads)
forward_ch.into { forward_ch_fasta ; forward_ch_unicycler ; forward_ch_masurca }

params.reverse_reads = "${launchDir}/DATA/atcc2001_S12_L001_R2_001.fastq.gz"
reverse_ch = Channel.fromPath(params.reverse_reads)
reverse_ch.into { reverse_ch_fasta ; reverse_ch_unicycler ; reverse_ch_masurca }

params.long_reads = "${launchDir}/DATA/atcc2001_nanopore.fastq.gz"
nanopore_ch = Channel.fromPath(params.long_reads)
nanopore_ch.into { nanopore_ch_fasta ; nanopore_ch_unicycler }

params.reference = "${launchDir}/DATA/GCF_000002545.3_ASM254v2_protein.faa" 
reference_ch = Channel.fromPath(params.reference)
reference_ch.into { unicycler_reference_ch ; smartdenovo_reference_ch ; masurca_reference_ch }


// Processes

process toFasta {

    input:
        file forward from forward_ch_fasta
        file reverse from reverse_ch_fasta
        file nanopore from nanopore_ch_fasta

    output:
        file "forward.fasta" into forward_fasta
        file "reverse.fasta" into reverse_fasta
        file "nanopore.fasta" into nanopore_fasta

    """
    gunzip ${forward} > forward.fastq
    paste - - - - < forward.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > forward.fasta
    gunzip ${reverse} > reverse.fastq
    paste - - - - < reverse.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > reverse.fasta
    gunzip ${nanopore} > nanopore.fastq
    paste - - - - < nanopore.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > nanopore.fasta
    """
}

nanopore_fasta.into { nanopore_fasta_smartdenovo ; nanopore_fasta_masurca }

process unicycler {

    label 'porky'

    input:
    file forward from forward_ch_unicycler
    file reverse from reverse_ch_unicycler
    file nanopore from nanopore_ch_unicycler

    output:
    file "assembly.fasta" into unicycler_assembly_ch

    """
    unicycler -t ${task.cpus} -o ${launchDir}/UNICYCLYER -1 ${forward} -2 ${reverse} -l ${nanopore}
    """
}

process smartdenovo {

    label 'porky'

    input:
    file forward from forward_fasta
    file reverse from reverse_fasta
    file nanopore from nanopore_fasta_smartdenovo

    output:
    file "assembly.fasta" into smartdenovo_assembly_ch

    """
    mkdir SMARTDENOVO; cd SMARTDENOVO/
    ### De novo assembly from Nanopore data ...
    smartdenovo.pl -p "smartdenovo" -t ${task.cpus} -k 21 -c 1 ${nanopore} > smartdenovo.mak
    make -f smartdenovo.mak
    ### Polish with Illumina reads twice ...
    ## 1
    # Map ...
    minimap2 -ax sr -t ${task.cpus} smartdenovo.dmo.cns ${forward} ${reverse} > smartdenovo_minimap2_1.sam
    # .sam to .bam
    samtools view -S -b smartdenovo_minimap2_1.sam > smartdenovo_minimap2_1.bam
    samtools sort -o smartdenovo_minimap2_1.sorted.bam smartdenovo_minimap2_1.bam
    samtools index smartdenovo_minimap2_1.sorted.bam
    # Polish ...    
    cp smartdenovo.dmo.cns smartdenovo.fasta
    java -Xmx16G -jar /TOOLS/pilon-1.23.jar --genome smartdenovo.fasta --frags smartdenovo_minimap2_1.sorted.bam
    mv pilon.fasta pilon_1.fasta
    ## 2
    # Map ...
    minimap2 -ax sr -t ${task.cpus} pilon_1.fasta ${forward} ${reverse} > smartdenovo_minimap2_2.sam    
    # .sam to .bam
    samtools view -S -b smartdenovo_minimap2_2.sam > smartdenovo_minimap2_2.bam
    samtools sort -o smartdenovo_minimap2_2.sorted.bam smartdenovo_minimap2_2.bam
    samtools index smartdenovo_minimap2_2.sorted.bam
    # Polish again ...
    java -Xmx16G -jar /TOOLS/pilon-1.23.jar --genome pilon_1.fasta --frags smartdenovo_minimap2_2.sorted.bam
    mv pilon.fasta assembly.fasta
    """
}

process masurca {

    label 'porky'

    input:
    file forward from forward_ch_masurca
    file reverse from reverse_ch_masurca
    file nanopore from nanopore_fasta_masurca

    output:
    file "assembly.fasta" into masurca_assembly_ch

    """
    mkdir MASURCA ; cd MASURCA/
    # Get configuration script, edit
    cp /TOOLS/masurca/sr_config_example.txt config.txt
    sed -i "s|PE= pe 500 50  /FULL_PATH/frag_1.fastq  /FULL_PATH/frag_2.fastq|PE= pe 500 50 ${forward} ${reverse}|g" config.txt
    sed -i "s|JUMP= sh 3600 200 |#JUMP= sh 3600 200 |g" config.txt
    sed -i "s|#NANOPORE=/FULL_PATH/nanopore.fa|NANOPORE=${nanopore}|g" config.txt
    sed -i "s|EXTEND_JUMP_READS=0|#EXTEND_JUMP_READS=0|g" config.txt
    sed -i "s|FLYE_ASSEMBLY=0|FLYE_ASSEMBLY=1|g" config.txt
    # Set up masurca and run    
    masurca config.txt
    ./assemble.sh
    """
}

process annotateUnicycler {

    label 'porky'

    input:
    file assembly from unicycler_assembly_ch
    file reference from unicycler_reference_ch

    output:
    file "augustus.hints.gtf" into unicycler_annotation_ch

    """
    mkdir /WORKSPACE/ANNOTATE
    mkdir /WORKSPACE/ANNOTATE/UNICYCLER ; cd /WORKSPACE/ANNOTATE/UNICYCLER/
    ## Simplify reference '.faa' headers ...
    simplifyFastaHeaders.pl ${reference} protein_aa_ protein_aa.fa protein_aa_header.map
    
    ## Soft mask the draft genome ...
    BuildDatabase -engine ncbi -name assembly ${assembly}
    RepeatModeler -database assembly -engine ncbi -pa ${task.cpus}
    RepeatMasker -pa ${task.cpus} -e ncbi -s -xsmall -lib assembly-families.fa ${assembly}
    
    ## Simplify masked assembly '.fasta' headers ...
    simplifyFastaHeaders.pl assembly.fasta.masked assembly_ assembly.fa assembly_header.map
    
    ## Run BRAKER
    braker.pl --genome=/WORKSPACE/assembly.fa --prot_seq=/WORKSPACE/protein_aa.fa --prg=gth --trainFromGth --softmasking --species=assembly --verbosity=4 --workingdir=/WORKSPACE/ANNOTATE/UNICYCLER/ --nocleanup --cores=${task.cpus}
    """
}

process annotateSmartdenovo {

    label 'porky'

    input:
    file assembly from smartdenovo_assembly_ch
    file reference from smartdenovo_reference_ch

    output:
    file "augustus.hints.gtf" into smartdenovo_annotation_ch

    """
    mkdir /WORKSPACE/ANNOTATE
    mkdir /WORKSPACE/ANNOTATE/SMARTDENOVO ; cd /WORKSPACE/ANNOTATE/SMARTDENOVO/
    ## Simplify reference '.faa' headers ...
    simplifyFastaHeaders.pl ${reference} protein_aa_ protein_aa.fa protein_aa_header.map
    
    ## Soft mask the draft genome ...
    BuildDatabase -engine ncbi -name assembly ${assembly}
    RepeatModeler -database assembly -engine ncbi -pa ${task.cpus}
    RepeatMasker -pa ${task.cpus} -e ncbi -s -xsmall -lib assembly-families.fa ${assembly}
    
    ## Simplify masked assembly '.fasta' headers ...
    simplifyFastaHeaders.pl assembly.fasta.masked assembly_ assembly.fa assembly_header.map
    
    ## Run BRAKER
    braker.pl --genome=/WORKSPACE/assembly.fa --prot_seq=/WORKSPACE/protein_aa.fa --prg=gth --trainFromGth --softmasking --species=assembly --verbosity=4 --workingdir=/WORKSPACE/ANNOTATE/SMARTDENOVO/ --nocleanup --cores=${task.cpus}
    """
}

process annotateMasurca {

    label 'porky'

    input:
    file assembly from masurca_assembly_ch
    file reference from masurca_reference_ch

    output:
    file "augustus.hints.gtf" into masurca_annotation_ch

    """
    mkdir /WORKSPACE/ANNOTATE
    mkdir /WORKSPACE/ANNOTATE/MASURCA ; cd /WORKSPACE/ANNOTATE/MASURCA/
    ## Simplify reference '.faa' headers ...
    simplifyFastaHeaders.pl ${reference} protein_aa_ protein_aa.fa protein_aa_header.map
    
    ## Soft mask the draft genome ...
    BuildDatabase -engine ncbi -name assembly ${assembly}
    RepeatModeler -database assembly -engine ncbi -pa ${task.cpus}
    RepeatMasker -pa ${task.cpus} -e ncbi -s -xsmall -lib assembly-families.fa ${assembly}
    
    ## Simplify masked assembly '.fasta' headers ...
    simplifyFastaHeaders.pl assembly.fasta.masked assembly_ assembly.fa assembly_header.map
    
    ## Run BRAKER
    braker.pl --genome=/WORKSPACE/assembly.fa --prot_seq=/WORKSPACE/protein_aa.fa --prg=gth --trainFromGth --softmasking --species=assembly --verbosity=4 --workingdir=/WORKSPACE/ANNOTATE/MASURCA/ --nocleanup --cores=${task.cpus}
    """
}

