__author__ = "Isra Nurhassen"

# imports
import glob
import subprocess
import shlex
import os
import os.path


# -------------- SET GLOBAL VARIABLES-------------- #


PATH_FASTQ = "/path/to/fastq/"
PATH_FASTQ_TR = "data/fastq/"
PATH_MAPPED = "data/mapped/"
PATH_SORTEDBAM = "data/sortedBAM/"
PATH_CALLS = "data/calls/"
PATH_LOGS = "data/logs/"
PATH_BED = "/path/to/bedfile.bed"
THREADS = 24
ID,= glob_wildcards(PATH_FASTQ + "{id}_R1_001.fastq.gz")


# --------------------- RULES ----------------------- #


rule all:
    input:
        "data/snakemake_job_complete.txt"


#  1 adapter and quality trimming
rule aq_trim:
    input:
        r1 = PATH_FASTQ + "{id}_R1_001.fastq.gz",
        r2 = PATH_FASTQ + "{id}_R2_001.fastq.gz"
    output:
        fastq_tr_R1 = temp(PATH_FASTQ_TR + "{id}_trimmed_R1.fastq"),
        fastq_tr_R2 = temp(PATH_FASTQ_TR + "{id}_trimmed_R2.fastq")
    params:
        t = "16"
    benchmark:
        "benchmarks/{id}/{id}.aq_trim.benchmark.txt"
    run:
        shell("fastp -i {input.r1} -o {output.fastq_tr_R1} -I {input.r2} -O {output.fastq_tr_R2} -w {params.t}")


#. 2 mapping reads against human reference genome
rule map_reads:
    input:
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        r1_tr = PATH_FASTQ_TR + "{id}_trimmed_R1.fastq",
        r2_tr = PATH_FASTQ_TR + "{id}_trimmed_R2.fastq"
    output:
        mapped = temp(PATH_MAPPED + "{id}.bam")
    params:
        t = THREADS
    log:
        PATH_LOGS + "{id}/{id}_bam.log"
    benchmark:
        "benchmarks/{id}/{id}.map_reads.benchmark.txt"
    run:
        shell("bwa-mem2.sse41 mem -t {params.t} {input.reference} {input.r1_tr} {input.r2_tr} 2> {log} | samtools view -bS -@{params.t} - > {output.mapped} ")


#  3 Sort BAM for calling
rule sort_bam:
    input:
        bam = rules.map_reads.output.mapped
    output:
        sorted = PATH_SORTEDBAM + "{id}_sorted.bam"
    params:
        t = THREADS
    log:
        PATH_LOGS + "{id}/03_{id}_sortedBAM.log"
    benchmark:
        "benchmarks/{id}/{id}.sort_bam.benchmark.txt"
    run:
        shell("samtools sort {input.bam} -o {output.sorted} -@{params.t} 2>&1 | tee {log}")


#  4 Marking duplicate reads in sorted BAM files 
rule mark_duplicates:
    input:
        bam = rules.sort_bam.output.sorted
    output:
        mDUP = temp(PATH_SORTEDBAM + "mDUP_{id}.bam"),
        metrics = temp(PATH_SORTEDBAM + "{id}_metrics.txt")
    log:
        PATH_LOGS + "{id}/04_{id}_markDuplicates.log"
    benchmark:
        "benchmarks/{id}/{id}.mark_duplicates.benchmark.txt"
    run:
        shell("picard -Xmx10g MarkDuplicates -INPUT {input.bam} -OUTPUT {output.mDUP} -M {output.metrics} 2>&1 | tee {log}")
    

#. 5 Add read groups
rule add_readGroups:
    input:
        bam = rules.mark_duplicates.output.mDUP
    output:
        RGbam = temp(PATH_SORTEDBAM + "RG_mDUP_{id}.bam")
    log:
        PATH_LOGS + "{id}/05_{id}_addReadgroups.log"
    # threads: THREADS
    benchmark:
        "benchmarks/{id}/{id}.add_readGroups.benchmark.txt"
    run:
        shell("picard -Xmx10g AddOrReplaceReadGroups -I {input.bam} -O {output.RGbam} -RGLB lib1 -RGPL ILLUMINA -SORT_ORDER coordinate -RGPU unit1 -RGSM {wildcards.id} 2>&1 | tee {log}")


#. 6 Recalibrate base with BQSR
#. BQSR corrects the quality scores in BAM file
rule recalibrate_base:
    input:
        bam = rules.add_readGroups.output.RGbam,
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        known_sites_snps = "/path/to/known_sites/known_sites/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz",
        known_sites_indels = "/path/to/known_sites/1000G_phase1.indels.hg19.sites.vcf.gz"
    output:
        recal_table = temp(PATH_SORTEDBAM + "recal_data_{id}.table")
    log:
        PATH_LOGS + "{id}/06_{id}_recalibrateBase.log"
    # threads: THREADS
    benchmark:
        "benchmarks/{id}/{id}.recalibrate_base.benchmark.txt"
    run:
        shell("gatk --java-options -Xmx4G BaseRecalibrator -I {input.bam} -R {input.reference} --known-sites {input.known_sites_snps} --known-sites {input.known_sites_indels} -O {output.recal_table} 2>&1 | tee {log}")


#. 7
rule applyBQSR:
    input:
        bam = PATH_SORTEDBAM+"RG_mDUP_{id}.bam",
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        recal_table = rules.recalibrate_base.output.recal_table
    output:
        bqsrBAM = PATH_SORTEDBAM + "BQSR_RG_mDUP_{id}.bam"
    log:
        PATH_LOGS + "{id}/07_{id}_applyBSQR.log"
    # threads: THREADS
    benchmark:
        "benchmarks/{id}/{id}.applyBQSR.benchmark.txt"
    run:
        shell("gatk --java-options -Xmx8G ApplyBQSR -R {input.reference} -I {input.bam} --bqsr-recal-file {input.recal_table} -O {output.bqsrBAM} 2>&1 | tee {log}")


#. 8 Variant calling with HaplotypeCaller
#  HaplotypeCaller in GVCF mode
rule call:
    input:
        bam = PATH_SORTEDBAM + "BQSR_RG_mDUP_{id}.bam",#rules.applyBQSR.output.bqsrBAM,
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        intervals = PATH_BED
    output:
        gvcf = "data/singlePatientsGVCF_fullGenes/{id}.g.vcf.gz"
    log:
        PATH_LOGS + "{id}/08_{id}_vcf.log"
    #threads: THREADS
    benchmark:
        "benchmarks/{id}/{id}.call.benchmark.txt"
    run:
        shell("gatk --java-options -Xmx8G HaplotypeCaller -R {input.reference} -I {input.bam} -O {output.gvcf} -ERC GVCF 2>&1 | tee {log}")
        #shell("gatk --java-options -Xmx8G HaplotypeCaller -R {input.reference} -I {input.bam} -O {output.gvcf} -L {input.intervals} -ERC GVCF 2>&1 | tee {log}"),


rule collection:
    input:
        dir = expand("data/singlePatientsGVCF_fullGenes/{id}.g.vcf.gz", id=ID)
    output:
         finished = temp("finished.txt")
    run:
        shell("touch {output.finished}")


#. 9
rule combineGVCFs:
    input:
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        finished = rules.collection.output.finished
    output:
        cohortGVCF = "data/combinedGVCF_fullGenes/cohort.g.vcf.gz"
    log:
         PATH_LOGS + "09_combineGVCF_vcf.log"
    #threads: THREADS
    benchmark:
        "benchmarks/combineGVCFs.benchmark.txt"
    run:
        allGVCF = ""
        for gvcf_file in ID:
            allGVCF = allGVCF + "-V " + "data/singlePatientsGVCF_fullGenes/" + gvcf_file + ".g.vcf.gz "
        shell("gatk --java-options -Xmx16g CombineGVCFs -R {input.reference} " + allGVCF + " -O {output.cohortGVCF} 2>&1 | tee {log}")

# 10
rule genotype:
    input:
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        cohortGVCF = rules.combineGVCFs.output.cohortGVCF
    output:
        vcf = PATH_CALLS + "allsamples.vcf.gz"
    log:
         PATH_LOGS + "2_genotype_vcf.log"
    #threads: THREADS
    benchmark:
        "benchmarks/genotype.benchmark.txt"
    run:
        shell("gatk --java-options -Xmx10g GenotypeGVCFs -R {input.reference} -V {input.cohortGVCF} -O {output.vcf} --annotation MappingQualityZero 2>&1 | tee {log}")



# 11 Seperate SNPs from Indels
rule selectVariants:
    input:
        vcf = rules.genotype.output.vcf
    output:
        snps = PATH_CALLS + "snps.vcf",
        indels = PATH_CALLS + "indels.vcf"
    log:
        snp = PATH_LOGS + "3_selectVariants_snp.log",
        indels = PATH_LOGS + "3_selectVariants_indels.log"
    #threads: THREADS
    benchmark:
        "benchmarks/selectVariants.benchmark.txt"
    run:
        shell("gatk --java-options -Xmx8G SelectVariants -V {input.vcf} -select-type SNP -O {output.snps} 2>&1 | tee {log.snp}"),
        shell("gatk --java-options -Xmx8G SelectVariants -V {input.vcf} -select-type INDEL -O {output.indels} 2>&1 | tee {log.indels}")



rule filterVariants:
    input:
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        snps_vcf = PATH_CALLS + "snps.vcf",
        indels_vcf = PATH_CALLS + "indels.vcf"
    output:
        snps_output = PATH_CALLS + "filtered_snp_output.vcf.gz",
        indels_output = PATH_CALLS + "filtered_indel_output.vcf.gz",
        dum_output = temp("dum_tamp.txt")
    log:
        PATH_LOGS + "13_applyVQSR_snp.log"
    #threads: THREADS
    benchmark:
        "benchmarks/filterVariants.benchmark.txt"
    run:
        shell("""gatk --java-options -Xmx8G VariantFiltration -R {input.reference} -V {input.snps_vcf} -O {output.snps_output} \
	-window 35 -cluster 3 \
	--filter-name "MQ" -filter "MQ < 40.0" \
	--filter-name "SOR3" -filter "SOR > 3.0" --filter-name "FS_60" -filter "FS > 60.0" \
	--filter-name "QUAL" -filter "QUAL < 30.0" \
	--filter-name "MQRankSum-12.5" -filter "MQRankSum < -12.5" \
	--filter-name "ReadPosRankSum-8" -filter "ReadPosRankSum < -8.0" \
	--filter-name "QD_2" -filter "QD < 2.0" \
	2>&1 | tee {log}"""),
        shell("""gatk --java-options -Xmx8G VariantFiltration -R {input.reference} -V {input.indels_vcf} -O {output.indels_output} \
        -window 35 -cluster 3 --filter-name "FS200" -filter "FS > 200.0" \
        --filter-name "QD2" -filter "QD < 2.0" \
        --filter-name "QUAL30" -filter "QUAL < 30.0" \
        --filter-name "ReadPosRankSum-20" -filter "ReadPosRankSum < -20.0" \
        --filter-name "SOR_gt_10" -filter "SOR > 10.0" \
         2>&1 | tee {log} """),
        shell("touch {output.dum_output}")

# 12
rule variantRecalibratorSNP:
    input:
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        snps_vcf = rules.filterVariants.output.snps_output,
        hapmap = "/path/to/hapmap/hapmap_3.3.hg19.sites.vcf.gz",
        omni = "/path/to/omni/1000G_omni2.5.hg19.sites.vcf.gz",
        g1k = "/path/to/g1k/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz",
        dbsnp = "/path/to/dbsnp/dbsnp_138.hg19.vcf.gz"
    output:
        vqsr_recal = "data/vqsr/output_vqsr_snps.recal",
        tranches = "data/vqsr/SNP.tranches"
    log:
        PATH_LOGS + "12_variantRecalibratorSNP.log"
    #threads: THREADS
    benchmark:
        "benchmarks/variantRecalibratorSNP.benchmark.txt"
    run:
        shell("""gatk --java-options -Xmx8G VariantRecalibrator -R {input.reference} -V {input.snps_vcf} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.g1k} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
        --max-gaussians 4 \
        -O {output.vqsr_recal} \
        --tranches-file {output.tranches} 2>&1 | tee {log}""")

# 13
rule applyVQSR_snp:
    input:
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        snps_vcf = rules.filterVariants.output.snps_output,
        recal = rules.variantRecalibratorSNP.output.vqsr_recal,
        tranches = rules.variantRecalibratorSNP.output.tranches
    output:
        snps_output = PATH_CALLS + "vqsr_snp_output.vcf.gz",
        dum_output = temp("dum_tamp.txt")
    log:
        PATH_LOGS + "13_applyVQSR_snp.log"
    #threads: THREADS
    benchmark:
        "benchmarks/applyVQSR_snp.benchmark.txt"
    run:
        shell("""gatk --java-options -Xmx8G ApplyVQSR -R {input.reference} -V {input.snps_vcf} -O {output.snps_output} \
        --truth-sensitivity-filter-level 99.7 \
        --tranches-file {input.tranches} \
        --recal-file {input.recal} \
        -mode SNP 2>&1 | tee {log} 
        """),
        shell("touch {output.dum_output}")

#14
rule variantRecalibratorINDELS:
    input:
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        indels_vcf = rules.filterVariants.output.indels_output,
        mills = "/path/to/mills/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
        dbsnp = "/path/to/dbsnp/dbsnp_138.hg19.vcf.gz",
        dum = rules.applyVQSR_snp.output.dum_output
    output:
        vqsr_recal = "data/vqsr_fullGenes/output_vqsr_indel.recal",
        tranches = "data/vqsr_fullGenes/INDEL.tranches"
    log:
        PATH_LOGS + "14_variantRecalibratorINDELS.log"
    #threads: THREADS
    benchmark:
        "benchmarks/variantRecalibratorINDELS.benchmark.txt"
    run:
        shell("""gatk --java-options -Xmx8G VariantRecalibrator -R {input.reference} -V {input.indels_vcf} \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 {input.mills} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an BaseQRankSum -an FS \
        -mode INDEL \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
        --max-gaussians 4 \
        -O {output.vqsr_recal} \
        --tranches-file {output.tranches} 2>&1 | tee {log} 
        """)


#. 15
rule applyVQSR_indel:
    input:
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        indels_vcf = rules.filterVariants.output.indels_output,
        recal = rules.variantRecalibratorINDELS.output.vqsr_recal,
        tranches = rules.variantRecalibratorINDELS.output.tranches
    output:
        indel_output = PATH_CALLS + "vqsr_indel_output.vcf.gz"
    log:
        PATH_LOGS + "15_applyVQSR_indel.log"
    #threads: THREADS
    benchmark:
        "benchmarks/applyVQSR_indel.benchmark.txt"
    run:
        shell("""gatk --java-options -Xmx8G ApplyVQSR -R {input.reference} -V {input.indels_vcf} -O {output.indel_output} \
        --truth-sensitivity-filter-level 95 \
        --tranches-file {input.tranches} \
        --recal-file {input.recal} \
        --create-output-variant-index true \
        -mode INDEL 2>&1 | tee {log}""")

#. 16 Filter genotype quality (GQ) and depth 
rule refinement:
    input:
        reference = "/mnt/db/genome/human/hg19_noalt_masked.fa",
        indels_vcf = rules.applyVQSR_indel.output.indel_output,
        snps_vcf = rules.applyVQSR_snp.output.snps_output
    output:
        snps_output = PATH_CALLS + "refined_vqsr_snp_output.vcf.gz",
        indels_output = PATH_CALLS +"refined_vqsr_indel_output.vcf.gz"
    log:
        PATH_LOGS + "16_refinement_variants.log"
    #threads: THREADS
    benchmark:
        "benchmarks/refineVariants.benchmark.txt"
    run:
        shell("""gatk --java-options -Xmx8G VariantFiltration -R {input.reference} -V {input.snps_vcf} -O {output.snps_output} \
        --genotype-filter-expression "GQ < 20" --genotype-filter-name "lowGQ" \
        --genotype-filter-expression "DP < 10" --genotype-filter-name "lowDP" \
         2>&1 | tee {log} 
        """),
        shell("""gatk --java-options -Xmx8G VariantFiltration -R {input.reference} -V {input.indels_vcf} -O {output.indels_output} \
       --genotype-filter-expression "GQ < 20" --genotype-filter-name "lowGQ" \
       --genotype-filter-expression "DP < 10" --genotype-filter-name "lowDP" \
         2>&1 | tee {log} 
        """)

rule end:
    input:
        rules.refinement.output
    output:
        "data/snakemake_job_complete.txt"
    shell:
        "touch {output}"


