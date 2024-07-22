PhantomBuster is a bioinformatical tool that removes phantom barcode combinations that occur due to single-nucleotide sequencing errors and index hopping.
It is written for lineage-tracing experiments and CRISPR-screens, but can be used for any experimental setups in which only barcodes and no genetic DNA is measured.

# Installation

PhantomBuster is available via [pypi](https://pypi.org/project/phantombuster/) and can be installed with standard python tools like pip or pipx.

    pipx install phantombuster

# QuickStart

PhantomBuster is a command line tool which can be run via the `phantombuster` command.
It consists of four main steps: (1) demultplexing, (2) error correction of random barcodes, (3) hopping removal and (4) thresholding.
For CRISPR-screens a separate script can calculate p-values for guides.

## Demultiplexing

PhantomBuster demultiplexes BAM or FASTQ files, extracts all specified barcodes while error correcting barcodes with known reference sequences.
For demultiplexing additional worker processes must be started.

```
phantombuster demultiplex [INPUTFILE] --outdir [DIR] --barcode-hierarchy-file [FILE] --regex-file [FILE] 
phantombuster worker --outdir [DIR]
```

INPUTFILE must be a csv file that lists all BAM and FASTQ files that are processed.

Example `input_files.csv`:
```
file
101.bam
```

The barcode hierarchy file is a csv file that lists all barcodes to be extracted.
The order of the barcodes creates a hierarchy, in which barcodes higher up the hierarchy are more general, while barcodes lower in the hierarchy are more specific.
The hierarchy is used in the second step error correction, in which two random barcode sequences are only compared, if all barcode sequences higher up the hierarchy are the same.
Example `barcode_hierarchy.csv`

```
barcode,type,referencefile,threshold,min_length,max_length
sample,reference,sample_barcodes.csv,auto,-,-
lib,reference,library_barcodes.csv,1,-,-
lid,random,-,-,50,50
```

The regex file is a csv file that specifies for each read region how to extract barcodes by a regular expression.

Example `regexes.csv`:
```
tag,regex
b2,"^[ACGTN]{3}(?P<sample>[ACGTN]{5})"
query,"(?P<lid>[ACGTN]{5,6}(?P<lib>ACGT|GTAC){s<=1}[ACGTN]+)"
```

The outdir is a directory that contains all output and temporary files.
The same out directory must be passed to all stages of phantombuster.

## Error Correction

The error correction step employs the UMI-tools error correction algorithm to error correct random barcode sequences.
For error correction additional worker processes must be started.

```
phantombuster error-correct --outdir [DIR] --barcode-hierarchy-file [FILE]
phantombuster worker --outdir [DIR]
```

The out dir and barcode hierarchy file must be the same as in the demultiplexing step.

## Index Hopping Removal

The index hopping removal step removes barcode combinations that likely arised due to index hopping.
For index hopping removal no worker processes need to be started.

```
phantombuster hopping-removal --outdir [DIR] [HOPPING_BARCODES]
```

The out dir must be the same as in the previous steps.
The barcodes to test must correspond to one or a combination of barcodes of the barcode hierarchy (`sample`).
Combinations can be given seperated by commas (`sample,lid`).
Multiple barcodes or combinations can be given and are then processed one after another (`sample,lid lib`)
It is recommended to test the combination of barcodes on the i5 index and then the combination of barcodes on the i7 index.

## Thresholding

Thresholding removes all barcode combinations with a read count below a user defined threshold.
Seperate Thresholds can be chosen for different values of barcodes, for example for each sample.
For tresholding no worker processes need to be started.

```
phantombuster threshold --outdir [DIR] --threshold-file [FILE]
```

Example `thresholds.csv`:
```
sample,threshold
sample1,10
```
