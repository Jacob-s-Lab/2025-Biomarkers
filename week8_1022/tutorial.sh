#!/usr/bin/sh
#SBATCH -A ACD114093          # Project name: The project number for our class.
#SBATCH -J fastqc             # Job name: You can change this to any name you prefer.
#SBATCH -p ngscourse          # Partition Name
#SBATCH -c 2                  # Number of CPU cores
#SBATCH --mem=13g             # Memory allocation
#SBATCH -o tutorial.out.log   # -o: This exports the out.log file, which will record the steps executed by the program.
#SBATCH -e tutorial.err.log   # -e: This exports the err.log file, which will record any failures; if not specified otherwise, both log files will be located in the current directory of the SH file.
#SBATCH --mail-user=          # Here, you can input your email. An email will be sent to you if the next line's conditions are met.
#SBATCH --mail-type=END       # Email will be sent when the job ends.

set -v -x
echo "start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"


# Please enter the R1 & R2 file name and your username
user=evelyn92
R1=SRR13076392_R1_5000
R2=SRR13076392_R2_5000
fastqdir=/work/${user}/result/fastqc
sample=SRR13076392_5000
 
sampleR1=${fastqdir}/${R1}.fastq
sampleR2=${fastqdir}/${R2}.fastq
echo "sampleR1 directory: " ${sampleR1}
echo "sampleR2 directory: " ${sampleR2}

## Create a new directory for FastQC report
tutorial_result=/work/${user}/result/tutorial_result
mkdir -p ${tutorial_result}
cd ${tutorial_result}
DIR_FQ=${tutorial_result}/tutorial_FQ
mkdir -p ${DIR_FQ}

echo "pwd for fastqc: " ${DIR_FQ}

## Set up the environment for running fastqc
module load biology
module load FastQC


## Analyzing your sample's sequence QC by fastqc
echo "fastqc start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

fastqc ${sampleR1} ${sampleR2} -o ${DIR_FQ}

echo "fastqc finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"


#------------fastqc finished----------------------


# ------------------------------------ #
# Please don't change the script below #
# ------------------------------------ #
## Reference: Homo_sapiens_assembly38.fasta
ref=/opt/ohpc/Taiwania3/pkg/biology/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta

echo "pwd for analysis: " ${tutorial_result}

echo "Analysis started"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

## Create a new directory for alignment
DIR_ALN=${tutorial_result}/tutorial_ALN
mkdir -p ${DIR_ALN}
cd ${DIR_ALN}

echo "pwd for alignment: "
pwd

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

bwa mem -M -R "@RG\tID:GP_${sample}\tSM:SM_${sample}\tPL:ILLUMINA" -t 40 -K 1000000 ${ref} ${sampleR1} ${sampleR2} > ${sample}.sam

echo "Mapping Reads: Finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"



######################################################
# Preparing for bam file  (sorting & indexing) #
######################################################
echo "preparing for bam file: start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

samtools view -@ 2 -S -b ${sample}.sam > ${sample}.bam
samtools sort -@ 2 ${sample}.bam -o ${sample}.sorted.bam
samtools index -@ 20 ${sample}.sorted.bam
echo "bam file has already prepared"
echo "$(date '+%Y-%m-%d %H:%M:%S')"




###################
# Mark duplicates #
###################
echo "Mark duplicates: Start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

PICARD=/work/opt/ohpc/Taiwania3/pkg/biology/Picard/picard_v2.27.4/share/picard-2.27.4-0/picard.jar
java -jar ${PICARD} MarkDuplicates -I ${sample}.sorted.bam -O ${sample}.sorted.markdup.bam -M ${sample}_markdup_metrics.txt --CREATE_INDEX true

echo "Mark duplicates: Finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"


## Create a new directory for variant calling
DIR_VC=${tutorial_result}/tutorial_VC
mkdir -p ${DIR_VC}
cd ${DIR_VC}

echo "pwd for variant calling: "
pwd

############################################
# Calling variants by GATK HaplotypeCaller #
############################################
echo "variant calling: start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

# set up the environment for variant calling
# GATK_PATH=/opt/ohpc/Taiwania3/pkg/biology/GATK/gatk_v4.2.3.0/gatk
module load Python
module load GATK/4.2.0.0

# python3 $(which gatk) HaplotypeCaller \
gatk HaplotypeCaller \
  -R /opt/ohpc/Taiwania3/pkg/biology/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta \
  -I ${DIR_ALN}/${sample}.sorted.markdup.bam \
  -O ${sample}.sorted.markdup.hc.vcf.gz \
  -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9  \
  -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17  \
  -L chr18 -L chr19 -L chr20 -L chr21 -L chr22

echo "variant calling: Finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"


## Create a new directory for annotation
DIR_ANNO=${tutorial_result}/tutorial_ANNO
mkdir -p ${DIR_ANNO}
cd ${DIR_ANNO}

echo "pwd for annotation: "
pwd

## Set up the environment and path for running VEP
VEP_PATH=/opt/ohpc/Taiwania3/pkg/biology/Ensembl-VEP/ensembl-vep/vep
VEP_CACHE_DIR=/opt/ohpc/Taiwania3/pkg/biology/DATABASE/VEP/Cache
VEP_FASTA=/opt/ohpc/Taiwania3/pkg/biology/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
BCFTOOLS=/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools

module load biology
module load Perl/5.28.1
export PATH=${PATH}:/opt/ohpc/Taiwania3/pkg/biology/HTSLIB/htslib_v1.13/bin:/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/samtools_v1.15.1/bin
set -euo pipefail


# split multiallelic
echo "split multiallelic: start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

${BCFTOOLS} norm -m -any ${DIR_VC}/${sample}.sorted.markdup.hc.vcf.gz \
    -Oz \
    -o ${sample}.sorted.markdup.hc.normed.vcf.gz
${BCFTOOLS} index -t -f ${sample}.sorted.markdup.hc.normed.vcf.gz

echo "Split multiallelic: Finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"


INPUT_VCF=${sample}.sorted.markdup.hc.normed.vcf.gz
SAMPLE_ID=${sample}.sorted.markdup.hc.normed.vcf.vep
echo "INPUT VCF directory: " ${INPUT_VCF}
echo "sample ID: " ${SAMPLE_ID}

#############################
# Variant annotation by VEP #
#############################
echo "VEP annotaion: start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

${VEP_PATH} --cache --offline \
    --cache_version 108 --dir_cache ${VEP_CACHE_DIR} \
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


echo "Finished all analysis"
echo "$(date '+%Y-%m-%d %H:%M:%S')"
