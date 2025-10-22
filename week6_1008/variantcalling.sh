#!/usr/bin/sh
#SBATCH -A ACD114093                     # Account name/project number
#SBATCH -J variantcalling                # Job name
#SBATCH -p ngscourse                     # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2                             # 使用的core數 請參考Queue資源設定
#SBATCH --mem=13g                        # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_vc.log                    # Path to the standard output file
#SBATCH -e err_vc.log                    # Path to the standard error output file
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


######################################################
# Preparing for bam file  (sorting & indexing) #
######################################################
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
# GATK_PATH=/opt/ohpc/Taiwania3/pkg/biology/GATK/gatk_v4.2.3.0
module load Biology
module load Python
module load GATK/4.2.3.0

gatk HaplotypeCaller \
	-R ${ref} \
	-I ${path2}/${sample}.sorted.markdup.bam \
	-O ${path3}/${sample}.HC.vcf.gz
echo "Variants calling: Finished"
echo "$(date '+%Y-%m-%d %H:%M:%S')"
