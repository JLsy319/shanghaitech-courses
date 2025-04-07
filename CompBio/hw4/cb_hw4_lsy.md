# Computational Biology Homework 4

## Q1

Illumina sequencing is high-throughput and ideal for whole genome sequencing, while Sanger sequencing is low-throughput, highly accurate, and best suited for targeted sequencing.

They share:

- Both rely on fluorescence to detect nucleotides. Sanger sequencing uses fluorescently labeled ddNTPs, while Illumina sequencing uses fluorescently labeled reversible terminators.

- Both methods typically require PCR amplification of the target DNA prior to sequencing.

## Q2

Differences:

|                     Illumina Sequencing                      |                          Nanopore Sequencing                           |
| :----------------------------------------------------------: | :--------------------------------------------------------------------: |
|            200-500 bp fragments, 50-300 bp reads             |     10-50kb or longer fragments, up to hundreds of kilobases reads     |
| Adapters contain primer binding sites for cluster generation | Adapters contain motor proteins for threading DNA through the nanopore |
|               Required to amplify the library                |              PCR is unnecessary for long-read sequencing               |
|      Sequencing by synthesis with fluorescent detection      |           Real-time sequencing via changes in ionic current            |

They share:

- Both support multiplexing through the use of barcodes.
- Both require quantification and quality assessment of the final library.

## Q3

Nanopore sequencing

The ionic current changes are sensitive to the chemical properties of the DNA bases, including modifications like methylation and other unknown modifications.

Also, there is no need for pre-treatment to detect modifications in nanopore sequencing. It can detect modification natively, which is helpful to detect unknown modifications.
