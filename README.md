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
| Option  | Use                                    | Default                   |
|:--------|:---------------------------------------|:--------------------------|
| --in    | Specify the input file in fasta format | You must specify an input |
| --out   | Specify the ouput folder name          | pipe_out                  |
| --motif | Specify the motif file to use          | ./data/motifs.txt         |

### Output files:
* (named output)
	* prokka-annotations
		* One folder for each annotation format, each containing annotation files for each processed sequence.
	* control-sequences
		* individual
			* Individual fasta files for each extracted control sequence.
		* all_sequences.fasta contains all the control sequences in one fasta file.
	* mast
		* Contains the MAST reports in html, json, and plaintext formats.

### To do:
It needs it to be laid out such that the ileu-tRNA comes before the 12S-rRNA, there's no error checking for anything, basically. If these features don't occur in the proper order, the output will be huge, so it's pretty easy to tell when this has happened.

Currently, if it doesn't find the tRNA for isoleucine, it will grab every base from 0 to the end, and also every base after the 12s rRNA. This means the "control region" will be larger than the entire genome. This isn't a *huge* problem, as it's pretty easy to tell when it's not working properly, but I should still get around to fixing it.

I think something similar happens if it doesn't find the 12S rRNA, but I forget.