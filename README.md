# MosMitCRT
MOSquito MITochondrial Control Region Trimmer

This software takes mitochondrial genomes, annotates them with Prokka, and extracts the control sequence (all the bases between the Isoleucine tRNA and the 12S rRNA that don't code for anything). The control sequence is then matched against a reference database of known species' control sequence motifs.

### Dependencies:
To run properly, this pipeline needs:
 * [Python 3](https://www.python.org/downloads/)
 * [Prokka](https://github.com/tseemann/prokka)
 * [MAST](http://meme-suite.org/doc/mast.html) (from the [MEME Suite](http://meme-suite.org/index.html))
 * [Beautiful Soup](https://www.crummy.com/software/BeautifulSoup/)
 * [lxml](https://lxml.de)
 * [bcbiogff](https://github.com/chapmanb/bcbb/tree/master/gff)

Recommended versions of theese dependencies are in the ```environment.yml``` file, for use with [Conda](https://docs.conda.io/en/latest/) virtual environments.

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
This command can also be run by using ```-profile test```.

No motif search, check reverses and rotated sequence handling:
```
nextflow run MosMitCRT --in MosMitCRT/testData/testv3.fasta --nomotif
```
This command can be run by using ```-profile test2```

### Versioning:

versioning: X.Y.Z
X: Major release version
  - Anything that adds an analysis step.
  - Things like adding a BLAST search, adding options for genome annotation programs, etc.
  - Changes the core functionality of the program.
Y: Minor release version
  - Changing the interface or adding smaller features, things that don't affect the core functionality of the program.
  - Giving a command line option to add Prokka options, changing output folder strucutre, etc.
  - Optionally different, but produces the same type of output and accepts the same inputs.
Z: Bugfix version
  - For when I make mistakes and have to push another version to change something that otherwise breaks the program.
  - Also used for back-end improvements that don't actually affect the usage or output.
