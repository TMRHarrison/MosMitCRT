#!/usr/bin/env nextflow

/* --in       file input(s)
 * --out      output folder name
 * --motif    notifs to scan for
 * --nomotif  skip motif search
 *
 */

def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run MosMitCRT --in data.fasta --out out_folder

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

// performs Prokka on the fasta files
process performProkka {
    // Publish to folders by extension, not by sequence.
    publishDir "${params.out}/prokka-annotations",
        mode: 'copy',
        saveAs: {fn ->
            ext = fn[fn.lastIndexOf(".")+1..-1]
            base = fn[fn.lastIndexOf("/")+1..(fn.lastIndexOf(".")-1)]
            "${ext}/${base}.${ext}"
        }

    input:
    file inp from sequences_ch // <--- --input

    output:
    file "**/*.gff" into annotatedSeqs_ch // ---> extractControlSeq
    file("**/*")

    """
    prokka --outdir prokka-out --force --prefix ${inp.baseName} --cpus 1 --gcode 5 --kingdom mitochondria ${inp}
    """


}

// extracts the control sequences
process extractControlSeq {
    publishDir "${params.out}/control-sequences/individual", mode: 'copy'

    input:
    file inp from annotatedSeqs_ch // <--- performProkka

    output:
    file "${inp.baseName}_controlSeq.fasta" into combineControlSeq_ch // ---> combineControlSeqs

    """
    extract-control.py --input ${inp} > ${inp.baseName}_controlSeq.fasta
    """

}

// concatenates the control sequences into a single file
process combineControlSeqs {
    publishDir "${params.out}/control-sequences", mode: 'copy'

    input:
    file "*.fasta" from combineControlSeq_ch.collect() // <--- extractControlSeqs

    output:
    file "all_sequences.fasta" into contSeqsCombined_ch // ---> split into (mast_ch -> findMotifs) and (annotateContReg_ch -> annotateMotifs)

    """
    cat *.fasta > all_sequences.fasta
    """

}

// when the control sequences file comes out, duplicate it
contSeqsCombined_ch
    .into {mast_ch; annotateContReg_ch}

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
