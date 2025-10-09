#!/usr/bin/sh
#SBATCH -A ACD114093           # Account name/project number
#SBATCH -J vep                 # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2                   # 使用的core數 請參考Queue資源設定
#SBATCH --mem=13g              # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_vep.log         # Path to the standard output file
#SBATCH -e err_vep.log         # Path to the standard error output file
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL,END


set -v -x
echo "start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"


# Please enter the R1 & R2 file name and your username
user=evelyn92
sampleR1=/work/${user}/result/fastqc/SRR13076392_S14_L002_R1_001.fastq.gz
sampleR2=/work/${user}/result/fastqc/SRR13076392_S14_L002_R2_001.fastq.gz
sample=SRR13076392
path1=/work/${user}/alignment/alignmentR
path2=/work/${user}/alignment/alignmentRM
path3=/work/${user}/variantcalling/variantcallingR


mkdir -p ${path1}
mkdir -p ${path2}
mkdir -p ${path3}

# ------------------------------------ #
# Please don't change the script below #
# ------------------------------------ #
## Reference: Homo_sapiens_assembly38.fasta
ref=/opt/ohpc/Taiwania3/pkg/biology/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
# Create the environment for alignment and variant calling
module load biology
module load BWA/0.7.17
module load SAMTOOLS/1.18 
set -euo pipefail




##############################
# Mapping reads with BWA-MEM #
##############################
echo "Mapping Reads: Start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

bwa mem -M -R "@RG\tID:GP_${sample}\tSM:SM_${sample}\tPL:ILLUMINA" -t 40 -K 1000000 ${ref} ${sampleR1} ${sampleR2} > ${path1}/${sample}.sam

echo "Mapping Reads: Finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"


################################################
# Preparing for bam file  (sorting & indexing) #
################################################
echo "preparing for bam file: start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

samtools view -@ 2 -S -b  ${path1}/${sample}.sam >  ${path1}/${sample}.bam
samtools sort -@ 2  ${path1}/${sample}.bam -o  ${path1}/${sample}.sorted.bam
samtools index -@ 20  ${path1}/${sample}.sorted.bam

echo "bam file has already prepared"
echo "$(date '+%Y-%m-%d %H:%M:%S')"


###################
# Mark duplicates #
###################
echo "Mark duplicates: Start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"
PICARD=/work/opt/ohpc/Taiwania3/pkg/biology/Picard/picard_v2.27.4/share/picard-2.27.4-0/picard.jar
java -jar ${PICARD} MarkDuplicates \
	-I  ${path1}/${sample}.sorted.bam \
	-O  ${path2}/${sample}.sorted.markdup.bam \
	-M  ${path2}/${sample}_markdup_metrics.txt \
	--CREATE_INDEX true
echo "Mark duplicates: Finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"


############################################
# Calling variants by GATK HaplotypeCaller #
############################################
echo "Variants calling: Start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

# set up the environment for variant calling
module load Python
module load GATK/4.2.0.0

gatk HaplotypeCaller \
	-R ${ref} \
	-I ${path2}/${sample}.sorted.markdup.bam \
	-O ${path3}/${sample}.HC.vcf.gz
echo "Variants calling: Finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

## Set up the environment and path for running VEP
VEP_PATH=/opt/ohpc/Taiwania3/pkg/biology/Ensembl-VEP/ensembl-vep/vep
VEP_CACHE_DIR=/opt/ohpc/Taiwania3/pkg/biology/DATABASE/VEP/Cache
VEP_FASTA=/opt/ohpc/Taiwania3/pkg/biology/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
BCFTOOLS=/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools

module load biology
module load Perl/5.28.1
module load Anaconda/Anaconda3
export PATH=${PATH}:/opt/ohpc/Taiwania3/pkg/biology/HTSLIB/htslib_v1.13/bin:/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/samtools_v1.15.1/bin
set -euo pipefail

# split multiallelic
echo "split multiallelic: start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"
${BCFTOOLS} norm -m -any ${sample}.HC.vcf.gz \
    -Oz \
    -o ${sample}.HC.normed.vcf.gz
${BCFTOOLS} index -t -f ${sample}.HC.normed.vcf.gz
echo "Split multiallelic: Finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"


INPUT_VCF=${sample}.HC.normed.vcf.gz
SAMPLE_ID=${sample}.HC.VEP
echo "INPUT VCF directory: " ${INPUT_VCF}
echo "sample ID: " ${SAMPLE_ID}
#############################
# Variant annotation by VEP #
#############################
echo "VEP annotaion: start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

${VEP_PATH} --cache --offline \
    --cache_version 108 \
    --dir_cache ${VEP_CACHE_DIR} \
    --assembly GRCh38 \
    --fasta ${VEP_FASTA} \
    --fork 4 \
    -i ${INPUT_VCF} \
       --check_existing \
    --af_gnomade \
    --af_gnomadg \
    --vcf \
    -o ${SAMPLE_ID}.vcf \
    --force_overwrite

echo "VEP annotation: Finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

# Generate tsv
echo -e "CHROM\tPOS\tREF\tALT\tDP\t$(${BCFTOOLS} +split-vep -l ${SAMPLE_ID}.vcf | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > ${SAMPLE_ID}.tsv
${BCFTOOLS} +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%CSQ\n' -A tab ${SAMPLE_ID}.vcf >> ${SAMPLE_ID}.tsv

awk 'NR==1 || ($1 ~ /^chr[1-9]$|^chr10$/) {
    gsub(/,.*$/, "", $6)
    print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $9 "\t" $10 "\t" $34 "\t" $44 "\t" $50 "\t" $54
}' ${SAMPLE_ID}.tsv > ${SAMPLE_ID}_filtered.tsv

echo "Format changing & column filtering: Finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

printf "##############################################################################\n"
printf "###              Work completed: $(date '+%Y-%m-%d %H:%M:%S')              ###\n"
printf "##############################################################################\n"
