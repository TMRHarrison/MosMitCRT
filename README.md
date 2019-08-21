# MosMitCRT
MOSquito MITochondrial Control Region Trimmer

This software takes mitochondrial genomes, annotates them with Prokka, and extracts the control sequence (all the bases between the Isoleucine tRNA and the 12S rRNA that don't code for anything). The control sequence is then matched against a reference database of known species' control sequence motifs.

### Dependencies:
To run properly, this pipeline needs:
 * [Python 3](https://www.python.org/downloads/)
 * [Prokka](https://github.com/tseemann/prokka)
 * [MAST](http://meme-suite.org/doc/mast.html) (from the [MEME Suite](http://meme-suite.org/index.html))

Prokka and MEME Suite can be easily installed via Conda:
```
conda install -c bioconda prokka
conda install -c bioconda meme
```

### Options:
| Option    | Use                                            | Default                                        |
|:----------|:-----------------------------------------------|:-----------------------------------------------|
| --in      | Specify the input file(s) in fasta format      | You must specify an input                      |
| --motif   | Specify the motif file to use                  | You must specify a motif file or use --nomotif |
| --out     | Specify the ouput folder name                  | pipe_out                                       |
| --nomotif | Override --motif and skip searching for motifs | Searches for motifs from the --motif file      |

### Output files:
* (named output)
	* prokka-annotations
		* One folder for each annotation format, each containing annotation files for each processed sequence.
	* control-sequences
		* individual
			* Individual fasta files for each extracted control sequence.
		* "all_sequences.fasta" contains all the control sequences in one fasta file.
	* mast
		* Contains the MAST reports in html, json, and plaintext formats.
	* "annotations.gff" contains annotations (and only annotations) in GFF3 format for the control sequences.
	* "annSeq.gff" contains annotations in GFF3 format, with a ##FASTA section containing all sequence data, as well.

### Test commad
This pipeline comes with some genomes from NCBI as test material, as well as some motifs to search for.

#### Test commands:
Normal usage:
```
nextflow run MosMitCRT --in "MosMitCRT/testData/mos*.fasta" --motif MosMitCRT/testData/motifs.txt
```

No motif search, check reverses and rotated sequence handling:
```
nextflow run MosMitCRT --in MosMitCRT/testData/testv3.fasta --nomotif
```

### To do:
Prokka doesn't analyze things as if they're circular, so if the sequence's 1 position is the middle of the 12S rRNA or the Ileu-tRNA, it won't recognise some or all of the feature. This could potentially be mitigated by running multiple
