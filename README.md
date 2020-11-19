# BAnalyzer - Barcode Analyzer for metabarcoding experiments

This pipeline provides some basic quality controls of a collection of barcode 
sequences for Metabarcoding experiments: taxid-wise dereplication, Hamming distance and sequence 
size distribution.

## Getting started

### Prerequisites 

BAnalyzer runs in a UNIX environment with BASH (tested on Debian GNU/Linux 10 
(buster)) and requires conda and an internet connection (at least for the first run).

### Installing

Start by getting a copy of this repository on your system, either by downloading and unpacking the archive, 
or using 'git clone':

```bash
cd path/to/repo/
git clone https://github.com/CVUA-RRW/BAnalyzer.git
```

Set up a conda environment containing snakemake, python and the pandas library and activate it:

```bash
conda create --name snakemake -c bioconda -c anaconda snakemake pandas
conda activate snakemake
```

### Sequence database

To run the pipeline you will need to provide a BLAST-formated reference sequence database.
If you already have a fasta file with your sequences follow the [BLAST documentation](https://www.ncbi.nlm.nih.gov/books/NBK279688/)
to know how to format it.

If you want to extract barcodes from a database of reference genomes you can check out 
our [RRW-PrimerBLAST](https://github.com/CVUA-RRW/RRW-PrimerBLAST) pipeline.

You will also need to provide the [taxdb](https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz) 
files available from the NCBI server.

### Running BAnalyzer

BAnalyzer should be run using the snakemake command-line application.
For this you will need to manually fill the config.yaml file with the paths to the required files.
You can also modify the parameters already present in the file.

Then run the pipeline with:

```bash 
snakemake -s /path/to/BAnalyzer/Snakefile --configfile path/to/config.yaml --use-conda --conda-prefix path/to/your/conda/envs
```

Consult [snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) for more details.

### Configuration file

The configuration file contains the following parameters:

```
# Fill in the path below with your own specifications:
workdir:                # Path to output directory
blast_db:               # Path to BLAST-formated database
taxdb:                  # Path to the folder containing the taxdb files

# Modify the parameters below:
trim_primers: False     # True to trim primers from sequences
primers: None           # Path to the fasta file containing primer sequences, required only if trim_primers is True
min_identity: 0.9       # Minimal identity level to compute Hamming distance (real between 0 and 1)
max_n: 3                # Maximum number of N nucleotides allowed per sequence, sequences with more will be discarded
```

If choosing to trim primer sequences from the barcodes, both the original and 
trimmed length will be shown in the report, but the Hamming-distance will be 
calculated only on the trimmed sequences.

### Taxid-wise clustering

Sequences within a Taxonomic node will be clustered prior to the alignement and calculation of the Hamming distance.
This allows to reduce redundancy in the database. Dereplication is performed by grouping 
sequences with 100% identity **within a taxonomic node**. 
Note that alignement of ambiguous nucleotides **never** incures a penalty. Therefore sequences
containing ambiguous nuclotides can be clustered with sequences that contain strict nucleotide (A, T, U, C, G). 

## Output

The pipeline produces an HTML report located in `workdir/reports` as well as different
CSV files that can be programmatically used for further analysis. 

## Credits

BAnalyzer is built with [Snakemake](https://snakemake.readthedocs.io/en/stable/) and uses:
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 
* [VSEARCH](https://github.com/torognes/vsearch)
* [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)

## Contributing

For new features or to report bugs please submit issues directly on the online repository.

## License

This project is licensed under a BSD 3-Clauses License, see the LICENSE file for details.

## Author

For questions about the pipeline, problems, suggestions or requests, feel free to contact:

Grégoire Denay, Chemisches- und Veterinär-Untersuchungsamt Rhein-Ruhr-Wupper 

<gregoire.denay@cvua-rrw.de>