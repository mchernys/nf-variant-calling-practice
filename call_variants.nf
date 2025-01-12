#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// bwa mem and gatk references, generated from the ensembl reference genome
params.reference_fasta = "/home/mchernys/references/ensembl/homo_sapiens/GRCh37/87/Homo_sapiens.GRCh37.dna.chromosome.7.fasta"
params.reference_dict = "/home/mchernys/references/ensembl/homo_sapiens/GRCh37/87/Homo_sapiens.GRCh37.dna.chromosome.7.dict"
params.reference_fai = "/home/mchernys/references/ensembl/homo_sapiens/GRCh37/87/Homo_sapiens.GRCh37.dna.chromosome.7.fasta.fai"
params.bwa_reference = tuple(params.reference_fasta, "${params.reference_fasta}.amb", "${params.reference_fasta}.ann", "${params.reference_fasta}.bwt", "${params.reference_fasta}.pac", "${params.reference_fasta}.sa")
params.gatk_reference = tuple(params.reference_fasta, params.reference_dict, "${params.reference_fasta}.fai")


// VCF files needed to build base quality recalibration model
// For GRCh37, they can be found at https://data.broadinstitute.org/snowman/hg19/variant_calling/vqsr_resources/Exome/v2/
// For GRCh38, they can be found at https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
params.broad_reference_dir="/home/mchernys/references/broad-institute/hg19/variant_calling/vqsr_resources/Exome/v2/"
params.broad_reference = tuple("${params.broad_reference_dir}/hapmap_3.3.b37.vcf.gz", "${params.broad_reference_dir}/hapmap_3.3.b37.vcf.gz.tbi",
                        "${params.broad_reference_dir}/1000G_omni2.5.b37.vcf.gz", "${params.broad_reference_dir}/1000G_omni2.5.b37.vcf.gz.tbi",
                        "${params.broad_reference_dir}/1000G_phase1.snps.high_confidence.b37.vcf.gz", "${params.broad_reference_dir}/1000G_phase1.snps.high_confidence.b37.vcf.gz.tbi",
                        "${params.broad_reference_dir}/dbsnp_138.b37.vcf.gz", "${params.broad_reference_dir}/dbsnp_138.b37.vcf.gz.tbi",
                        "${params.broad_reference_dir}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz", "${params.broad_reference_dir}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz.tbi")

// primers to trim
params.forward_primers = "data/forward_primers.fasta"
params.reverse_primers = "data/reverse_primers.fasta"
// contains folders with fastq files R1.fastq.gz and R2.fastq.gz
params.fastq_dir = 'data/raw'
// names of folders in fastq_dir to use. These names will be used as sample names
params.samples = "data/samples.txt"
// where to publish bam files
params.alignment_dir = "data/GRCh37/alignment/"
// where to publish vcf files
params.variant_dir = "data/GRCh37/variant_calling/"


params.cores = 10

process cutadapt {
    maxForks 2
    tag "$sample"

    input:
        tuple val(sample), path(r1_fastq), path(r2_fastq)
        path forward_primers
        path reverse_primers

    output:
        tuple val(sample), path("cutadapt-$sample-R1.fastq.gz"), path("cutadapt-$sample-R2.fastq.gz")

    script:
        """
        forward_primer_list=\$(sed -n -e '0~2p' $forward_primers | sed -z 's/\\n/ -g /g')
        reverse_primers_list=\$(sed -n -e '0~2p' $reverse_primers | sed -z 's/\\n/ -G /g')
        cutadapt -j $params.cores -g \${forward_primer_list} -G \${reverse_primers_list} -o cutadapt-$sample-R1.fastq.gz -p cutadapt-$sample-R2.fastq.gz $r1_fastq $r2_fastq
        """
}

process bwa_mem {
    maxForks 1

    tag "$sample"

    input:
        tuple val(sample), path("cutadapt-$sample-R1.fastq.gz"), path("cutadapt-$sample-R2.fastq.gz")
        tuple path(reference_fasta), path(amb), path(ann), path(bwt), path(ac), path(sa) // bwa index files

    output:
        tuple val(sample), path("${sample}.sam.gz")

    script:
        """
        bwa mem -M -t $params.cores $reference_fasta cutadapt-$sample-R1.fastq.gz cutadapt-$sample-R2.fastq.gz | gzip -3 > ${sample}.sam.gz
        """
}

process sam2bam {
    maxForks 4

    tag "$sample"

    input:
        tuple val(sample), path("${sample}.sam.gz")

    output:
        tuple val(sample), path("${sample}.bam")

    script:
        """
        samtools view -S -b ${sample}.sam.gz > ${sample}.bam
        """
}


process groupSortIndex {
    maxForks 4

    tag "$sample"

    input:
        tuple val(sample), path("${sample}.bam")

    output:
        tuple val(sample), path("${sample}.sorted.bam")

    script:
        """
        gatk AddOrReplaceReadGroups I=${sample}.bam O=${sample}.named.bam \
                RGLB=lib1 \
                RGPL=ILLUMINA \
                RGPU=unit1 \
                RGSM=${sample} \
                RGID=${sample}
        samtools sort ${sample}.named.bam -o ${sample}.sorted.bam
        samtools index ${sample}.sorted.bam
        """
}

process BQSR {
    maxForks 4
    publishDir params.alignment_dir

    tag "$sample"

    input:
        tuple val(sample), path("${sample}.sorted.bam")
        tuple path(reference_fasta), path(dict), path(fai) // gatk reference files
        tuple path(hapmap), path(hapmap_index), path(omni), path(omni_index), path(phase1), path(phase1_index), path(dbsnp), path(dbsnp_index), path(indels), path(indels_index) // VCF files used for BQSR/VQSR

    output:
        tuple val(sample), path("${sample}.BQSR.bam")

    script:
        """
        gatk --java-options "-Xmx4g" BaseRecalibrator \
            -I ${sample}.sorted.bam \
            -R ${reference_fasta} \
            -O ${sample}_BQSR_recal.table \
            --known-sites ${hapmap} \
            --known-sites ${omni} \
            --known-sites ${phase1} \
            --known-sites ${dbsnp} \
            --known-sites ${indels}

        gatk ApplyBQSR \
            -R ${reference_fasta} \
            -I ${sample}.sorted.bam \
            --bqsr-recal-file ${sample}_BQSR_recal.table \
            -O ${sample}.BQSR.bam
        """
}

process haplotypeCaller {
    maxForks 2

    tag "$sample"

    input:
        tuple val(sample), path("${sample}.BQSR.bam")
        tuple path(reference_fasta), path(dict), path(fai) // gatk reference files

    output:
        path "${sample}.g.vcf.gz"

    script:
        """
        gatk --java-options "-Xmx4g" HaplotypeCaller  \
            --output-mode EMIT_ALL_ACTIVE_SITES \
            --dont-use-soft-clipped-bases \
            --all-site-pls true \
            --native-pair-hmm-threads ${params.cores} \
            -R ${reference_fasta} \
            -I ${sample}.BQSR.bam \
            -O ${sample}.g.vcf.gz \
            -ERC GVCF;
        """
}

process genotypeGVCFs {

    publishDir params.variant_dir

    input:
        val gvcfs
        tuple path(reference_fasta), path(dict), path(fai) // gatk reference files

    output:
        tuple path("cohort.vcf.gz"), path("cohort.vcf.gz.tbi")

    script:
        """
        vcf_combine_arg=\$(echo $gvcfs | sed -e 's/, / --variant /g' -e 's/\\[//g' -e 's/\\]//g')
        echo \${vcf_combine_arg}
        gatk GenomicsDBImport -V \${vcf_combine_arg} --genomicsdb-workspace-path genomicsdb --intervals 7
        gatk GenotypeGVCFs --reference $reference_fasta --variant gendb://genomicsdb --output cohort.vcf.gz
        """
}

process VQSR {

    publishDir params.variant_dir

    input:
        tuple path("cohort.vcf.gz"), path("cohort.vcf.gz.tbi")
        tuple path(reference_fasta), path(dict), path(fai) // gatk reference files
        tuple path(hapmap), path(hapmap_index), path(omni), path(omni_index), path(phase1), path(phase1_index), path(dbsnp), path(dbsnp_index), path(indels), path(indels_index) // VCF files used for BQSR/VQSR

    output:
        tuple path("cohort.VQSR.vcf.gz"), path("cohort.VQSR.vcf.gz.tbi")

    script:
        """
        # calibrate for SNPs
        gatk VariantRecalibrator \
            -R ${reference_fasta} \
            -V cohort.vcf.gz \
            --max-gaussians 4 \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${phase1} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode SNP \
            -O cohort.SNP.recal \
            --tranches-file cohort.SNP.tranches \
            --rscript-file cohort.SNP.plots.R

        # calibrate for INDELs
        gatk VariantRecalibrator \
            -R ${reference_fasta} \
            -V cohort.vcf.gz \
            --max-gaussians 4 \
            --resource:mills,known=false,training=true,truth=true,prior=12.0 ${indels} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
            -mode INDEL \
            -O cohort.INDEL.recal \
            --tranches-file cohort.INDEL.tranches \
            --rscript-file cohort.INDEL.plots.R

        # apply SNP score cutoff to filter variants based on recalibration table (model)
        gatk ApplyVQSR \
            -R ${reference_fasta} \
            -V cohort.vcf.gz \
            -O SNP-cohort.vcf.gz \
            --exclude-filtered true \
            --truth-sensitivity-filter-level 100 \
            --tranches-file cohort.SNP.tranches \
            --recal-file cohort.SNP.recal \
            -mode SNP

        # apply INDEL score cutoff to filter variants based on recalibration table (model)
        gatk ApplyVQSR \
            -R ${reference_fasta} \
            -V SNP-cohort.vcf.gz \
            -O cohort.VQSR.vcf.gz \
            --exclude-filtered true \
            --truth-sensitivity-filter-level 100 \
            --tranches-file cohort.INDEL.tranches \
            --recal-file cohort.INDEL.recal \
            -mode INDEL
        """

}

workflow {
    println "Fastq files to process: "
    def samples_ch = channel.fromPath(params.samples).splitText(){it.strip()}.map(sample -> tuple(sample, "$params.fastq_dir/$sample/R1.fastq.gz", "$params.fastq_dir/$sample/R2.fastq.gz")).view()
    samples_ch = cutadapt(samples_ch, params.forward_primers, params.reverse_primers)
    samples_ch = bwa_mem(samples_ch, params.bwa_reference)
    samples_ch = sam2bam(samples_ch)
    samples_ch = groupSortIndex(samples_ch)
    samples_ch = BQSR(samples_ch, params.gatk_reference, params.broad_reference)
    samples_ch_list = haplotypeCaller(samples_ch, params.gatk_reference).collect().view()
    cohort_vcf = genotypeGVCFs(samples_ch_list, params.gatk_reference).collect()
    VQSR(cohort_vcf, params.gatk_reference, params.broad_reference)
}
