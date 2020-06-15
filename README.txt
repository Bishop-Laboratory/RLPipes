

Illumina iGenomes

The iGenomes are a collection of sequence and annotation files for commonly 
analyzed genomes. Each iGenome contains data for one species, downloaded from 
one source (UCSC, NCBI, or Ensembl), for one genomic build. 


*************
CHANGE POLICY
*************

Please see the CHANGE_LOG.txt file for all changes to the iGenomes.

The reference sequence of a build will not change. The only exception is if 
there has been an error in the sequencing, as reported by the source (eg NCBI).
In this case, the iGenome reference sequence will be corrected, and an entry
added to the CHANGE_LOG.txt file.

Annotation files of a build will be updated approximately every year.
The most current annotation files are in the Annotation directory. In this
directory, there is a README.txt file that contains the date of the most
recent update for that build. Also in the Annotation directory, the Archives
directory contains annotation data from prior updates. When annotation updates
occur, entries will be added to CHANGE_LOG.txt.

No files or directory will be modified or deleted without an accompanying
entry in the CHANGE_LOG.txt file. Some files may be added, such as new iGenome
species or builds, and these additions will be noted in CHANGE_LOG.txt.


************
INSTALLATION
************

The sequence and annotation files for each iGenome are compressed in 
a .tar.gz file. After downloading this file for a given genome, uncompress 
it with the following example unix/linux command (substituting 'iGenome' 
for the filename you downloaded):
    tar -zxvf iGenome.tar.gz

Each iGenome will uncompress with the directory structure described below.
Multiple iGenomes can be uncompressed in the same directory, creating
an integrated directory structure by Species/Source/Build. 
For example, uncompressing two different builds of the UCSC human genome:
    tar -zxvf Homo_sapiens_UCSC_hg18.tar.gz
    tar -zxvf Homo_sapiens_UCSC_hg19.tar.gz
will yield the following nested directory structure:
    Homo_sapiens/
    UCSC/
    hg18/ hg19/

To uncompress many iGenomes concurrently, use this command:
    for f in *.tar.gz; do tar -zxvf "$f"; done


****************
iGENOME CONTENTS
****************

Each iGenomes has the following nested directory structure:
    Species/
    Source/
    Build/
    Annotation/ Sequence/

Sequence/ contains subdirectories with genome sequences in various file
formats.

Sequence/AbundantSequences/ contains known abundant sequences for a genome. 
Every iGenome will contain phiX and sequencing adaptor sequences, in both fasta 
and squashed format. This directory is useful for RNA analysis.

Sequence/BowtieIndex/ contains an index of the whole genome for use with the
Bowtie aligner, which is also used by TopHat/Cufflinks. 

Sequence/Bowtie2Index/ contains an index of the whole genome for use with the
Bowtie2 aligner, which is also used by TopHat2.

Sequence/BWAIndex/ contains an index of the whole genome for use with the
BWA aligner. Since the BWA index scheme changed at version 0.6.0, both
versions 0.5.x and 0.6.0 (and higher) are included.

Sequence/Chromosomes/ contains individual chromosome sequences in fasta format.
Chromosomes are named according to their download source. Genome sequences from
UCSC are repeat-masked with lower-case characters; those from NCBI are not
repeat-masked; and those from Ensembl are repeat-masked with 'N' characters.
Masking was performed by the download source prior to download for the iGenomes.

Sequence/WholeGenomeFasta/ contains a multi-fasta files that contains all
chromosome sequences. This directory also contains samtools index file and
a picard sequence dictionary.


Annotation/ contains subdirectories with annotation files. Since not all
file types are present for each genome, only those files available from each 
source are provided.

Annotation/Genes/ contains files related to gene annotation if they are
available from a download source. 'Seq_gene.md.gz' files are downloaded from
NCBI, 'refFlat' files from UCSC, and 'GTF' files are from Ensembl. To ensure
that consistent formats are available, these files are converted such that all
genomes (with gene annotation) have refFlat.txt.gz and genes.gtf files. For
example, refFlat files downloaded from UCSC are also converted to gtf format.
Genes/ may also contain: ChromInfo.txt, refGene.txt, cytoBand.txt, kgXref.txt, 
knownGene.txt, knownToRefSeq.txt, refMrna.fa, refSeqSummary.txt,

A README.txt file in Annotation/ contains the date when the most recent 
annotation files were downloaded. For Ensembl genomes, this file also contains 
the Ensembl release number of the most recent annotation files. 

Annotation/SmallRNA/ contains species-specific files downloaded from miRBase,
for species where data is available.

Annotation/Variation/ contains files related to variation within the species.
Variation is currently supported for only a few species. The file formats are 
'chromosome_reports' (NCBI) and source-specific mysql table dumps in 
tab-delimited text format (UCSC).



**********************
ADDITIONAL INFORMATION
**********************

Symbolic links have been included within iGenome directories to reduce
total file sizes.

Several files have been modified so they will work with Tophat/Cufflinks. 
See the README.txt in the Annotation directory of each iGenome for details of 
any modifications.


