# ![qiubio-nf-core/atacseq](assets/atacseqlogo.png)

## Introduction

**qiubio/atacseq** is a bioinformatics analysis pipeline used for Assay for Transposase-Accessible Chromatin using sequencing (ATAC-seq) data based on [nfcore/atacseq](https://nf-co.re/atacseq).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

This pipeline will do qc analysis by ATACseqQC, generate the UCSC genome browser track hub and metagene analysis resuls in addition to original output of **nfcore/atacseq** pipeline.

## Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Alignment ([`BWA`](https://sourceforge.net/projects/bio-bwa/files/))
4. Mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
5. Merge alignments from multiple libraries of the same sample ([`picard`](https://broadinstitute.github.io/picard/))
    1. Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
    2. Filtering to remove:
        * reads mapping to mitochondrial DNA ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads mapping to blacklisted regions ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/))
        * reads that are marked as duplicates ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that arent marked as primary alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that are unmapped ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that map to multiple locations ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads containing > 4 mismatches ([`BAMTools`](https://github.com/pezmaster31/bamtools))
        * reads that are soft-clipped ([`BAMTools`](https://github.com/pezmaster31/bamtools))
        * reads that have an insert size > 2kb ([`BAMTools`](https://github.com/pezmaster31/bamtools); *paired-end only*)
        * reads that map to different chromosomes ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
        * reads that arent in FR orientation ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
        * reads where only one read of the pair fails the above criteria ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
    3. Alignment-level QC and estimation of library complexity ([`picard`](https://broadinstitute.github.io/picard/), [`Preseq`](http://smithlabresearch.org/software/preseq/))
    4. Create normalised bigWig files scaled to 1 million mapped reads ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
    5. Generate gene-body meta-profile from bigWig files ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html))
    6. Calculate genome-wide enrichment ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html))
    7. Call broad/narrow peaks ([`MACS2`](https://github.com/taoliu/MACS))
    8. Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
    9. Create consensus peakset across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
    10. Count reads in consensus peaks ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
    11. Differential accessibility analysis, PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
    12. Generate ATAC-seq specific QC html report ([`ataqv`](https://github.com/ParkerLab/ataqv))
6. Merge filtered alignments across replicates ([`picard`](https://broadinstitute.github.io/picard/))
    1. Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
    2. Remove duplicate reads ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    3. Create normalised bigWig files scaled to 1 million mapped reads ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
    4. Call broad/narrow peaks ([`MACS2`](https://github.com/taoliu/MACS))
    5. Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
    6. Create consensus peakset across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
    7. Count reads in consensus peaks relative to merged library-level alignments ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
    8. Differential accessibility analysis, PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
7. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation ([`IGV`](https://software.broadinstitute.org/software/igv/)).
8. Present QC for raw read, alignment, peak-calling and differential accessibility results ([`ataqv`](https://github.com/ParkerLab/ataqv), [`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

## Installation by conda

```bash
conda update conda
wget https://raw.githubusercontent.com/jianhong/atacseq/master/environment.yml
conda env create -n atacflow -f environment.yml
rm environment.yml
conda activate atacflow
nextflow pull jianhong/atacseq -r dev
srun --mem 60G -c 2 nextflow run jianhong/atacseq -profile test -r dev
```

## Update

```bash
conda activate atacflow
nextflow pull jianhong/atacseq
```

## Remove

```bash
conda activate atacflow
nextflow drop jianhong/atacseq
conda deactivate
conda remove --name atacflow --all
conda info --envs
```

## design table

| group | replicate | fastq_1 | fastq_2 | track_color | track_group |
|-------|-----------|---------|---------|-------------|-------------|
| WT | 1 | fastq/WT1.fastq.gz| | #E69F00 | SAMPLE |
| WT | 2 | fastq/WT2.fastq.gz| | #E69F00 | SAMPLE |
| KD | 1 | fastq/KD1.fastq.gz| | #0000FF | SAMPLE |
| KD | 2 | fastq/KD2.fastq.gz| | #0000FF | SAMPLE |
| Input | 1 | fastq/KD1.fastq.gz| | #000000 | SAMPLE |
| Input | 2 | fastq/KD2.fastq.gz| | #000000 | SAMPLE |

## metagene analysis

```
nextflow run jianhong/atacseq -profile test -resume --genomicElements beds/*.bed
```

## Get help

Please create an issue to submit your questions.



