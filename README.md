# nf-variant-calling-practice

This is a practice project to learn variant calling using Nextflow, done according the the GATK4 best practices. The pipeline is designed to be run on a cohort of human samples on GRCh37, but can be easily modified to run on a single sample.

The pipeline trims forward and reverse primers from the reads, aligns them using bwa, and calls variants using GATK4.

## Dependencies
- Quality control: [FastQC v0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC v1.14](https://multiqc.info/)
- Primer trimming: [cutadapt v3.2](https://cutadapt.readthedocs.io/en/stable/installation.html)
- Alignment: [bwa v0.7.17-r1188](https://github.com/lh3/bwa), [samtools v1.16.1](http://www.htslib.org/download/)
- Variant calling: [GATK v4.3.0.0](https://github.com/broadinstitute/gatk/releases)


## Building references
```
bowtie2-build Homo_sapiens.GRCh37.dna.chromosome.7.fa Homo_sapiens.GRCh37.dna.chromosome.7
samtools faidx Homo_sapiens.GRCh37.dna.chromosome.7.fa
gatk CreateSequenceDictionary -R Homo_sapiens.GRCh37.dna.chromosome.7.fa
```

## Quality control

- Utilize FastQC to create quality reports on raw fastq.gz reads and combine these reports using multiQC.

```
for d in ../data/raw/*/; do
   donor="$(basename "$d")"
   ln -sf $(realpath $d/R1.fastq.gz) ../data/fastqc/${donor}-R1.fastq.gz;
   ln -sf $(realpath $d/R2.fastq.gz) ../data/fastqc/${donor}-R2.fastq.gz;
done

fastqc --threads 10 ../data/fastqc/*
multiqc ../data/fastqc/ --outdir ../data/fastqc/
```

## Notes
- Variant calling performed according to [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode).
- The trickiest part of the pipeline is the ```gatk VariantRecalibrator``` step, which trains a model to flexibly filter variants. See [the VQSR docs](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-)
- Variant recalibration resources for GRCh37 can be found [on the broad institute website](https://data.broadinstitute.org/snowman/hg19/variant_calling/vqsr_resources/Exome/v2/)
- Variant recalibration resources for GRCh38 is linked [in the GATK documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)
- ```bwa mem``` aligns raw fastq files, do not merge them beforehand.
- We add library information to bam files in order to be able to combine data in the cohort vcf file.
- ```gatk HaplotypeCaller``` realigns reads within an "active" area where there is suspected variation, so the variants might not line up perfectly with the bam alignments.
- SNP variant recalibration model built using ```gatk VariantRecalibrator```
- INDEL variant recalibration model built using ```gatk VariantRecalibrator```
- Apply SNP and INDEL variant recalibration model filtering using ```gatk ApplyVQSR```
- use ```samtools faidx``` to select vcf region, then ```bcftools consensus``` to draw a consensus
