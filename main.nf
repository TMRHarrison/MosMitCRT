#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run TMRHarrison/MosMitCRT --in data.fasta --out out_folder

    Mandatory arguments:
      --in                      Path to input fasta file(s).
      --motif                   Motif file location.

    Other options:
      --out                     Output folder name. Defaults to "pipe_out"
      --nomotif                 Overrides --motif and skips MAST and subsequent steps.


    """.stripIndent()
}

params.in = null
params.out = "pipe_out"
params.motif = null
params.help = null
params.nomotif = null

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}


if (params.in == null) {
    log.info"You must specify an input FASTA file with \"--in <file>\"."
    exit 1
}
if (params.motif == null && !params.nomotif) {
    log.info"You must supply a motif file with \"--motif <file>\" or use the \"--nomotif\" flag to skip motif searching."
    exit 1
}


// Split the sequences into their own processes.
sequences_ch = Channel.fromPath(params.in)
                    .splitFasta(by: 1, file: true)
                    .dump()

// load up the reference motifs
if (!params.nomotif) {
    referenceMotifs_ch = Channel.fromPath(params.motif)
} else {
    referenceMotifs_ch = Channel.empty()
}

// appends the sequence ID to the filename. Only takes the ID up until the first non-alphanumeric character that isn't '.', '_', or '-'.
// This behaviour is to prevent illegal filenames.
// The rest of the filename is maintained to ensure the output won't have any overlapping names.
process giveFileNameFastaID {
    input:
    file inp from sequences_ch// <--- --in

    output:
    file "*_${inp}" into renamedSequences_ch // ---> performProkka

    """
    fastaID=\$(awk -F "[^a-zA-Z0-9\\._-]" '/^>/ {print \$2; exit}' ${inp})
    cp $inp \${fastaID}_${inp}
    """
}

// performs Prokka on the fasta files
process performProkka {
    // Publish to folders by extension, not by sequence.
    publishDir "${params.out}/prokka-annotations",
        mode: 'copy',
        saveAs: {filename ->
            fileOut = file(filename)
            "${fileOut.getExtension()}/${fileOut.getName()}"
        }

    input:
    file inp from renamedSequences_ch // <--- giveFileNameFastaID

    output:
    file "**/*.gff" into annotatedSeqs_ch // ---> extractControlSeq
    file("**/*")

    """
    prokka --outdir prokka-out --force --prefix ${inp.baseName} --cpus ${task.cpus} --gcode 5 --kingdom mitochondria ${inp}
    """


}

// extracts the control sequences
process extractControlSeq {
    publishDir "${params.out}/control-sequences/individual", mode: 'copy'

    input:
    file inp from annotatedSeqs_ch // <--- performProkka

    output:
    file fasta into combineControlSeq_ch // ---> combineControlSeqs

    script:
    fasta = "${inp.baseName}_controlSeq.fasta"

    """
    extract-control.py --input ${inp} > ${fasta}
    """

}

// concatenates the control sequences into a single file
process combineControlSeqs {
    publishDir "${params.out}/control-sequences", mode: 'copy'

    input:
    file "*.fasta" from combineControlSeq_ch.collect() // <--- extractControlSeqs

    output:
    file "all_sequences.fasta" into mast_ch, annotateContReg_ch // ---> split into (mast_ch -> findMotifs) and (annotateContReg_ch -> annotateMotifs)

    """
    cat *.fasta > all_sequences.fasta
    """

}

// Uses mast to find known motifs
process findMotifs {
    publishDir params.out, mode: 'copy'

    input:
    file inp from mast_ch // <--- combineControlSeqs
    file ref from referenceMotifs_ch // <--- --motifs (command line or defualt file)

    output:
    file "mast/mast.xml" into annotateMotifs_ch // ---> annotateMotifs
    file "mast/**" // ---> output

    """
    mast -o ./mast ${ref} ${inp}
    """

}

// turns the mast annotations into .gff format, and makes a second file that has annotations and the sequences.
process annotateMotifs {
    publishDir params.out, mode: 'copy'

    input:
    file inp from annotateMotifs_ch // <--- findMotifs
    file seqs from annotateContReg_ch // <--- input

    output:
    file "annotations.gff" // ---> output
    file "annSeq.gff" // ---> output

    """
    mast-annotate.py --mast ${inp} > annotations.gff && \
    bind-gff-to-fasta.py --gff annotations.gff --fasta ${seqs} > annSeq.gff
    """

}
