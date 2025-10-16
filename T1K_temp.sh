#!/usr/bin/sh
#SBATCH -A ACD114093 
#SBATCH -J person_name_T1K
#SBATCH -p ngscourse
#SBATCH -c 2
#SBATCH --mem=13g
#SBATCH -o out_person_name_T1K.log
#SBATCH -e err_person_name_T1K.log
#SBATCH --mail-user
#SBATCH --mail-type=FAIL

user=evelyn92
mkdir -p /work/${user}/HLA_typing
sampleR1=/work/${user}/result/fastqc/SRR13076392_S14_L002_R1_001.fastq.gz
sampleR2=/work/${user}/result/fastqc/SRR13076392_S14_L002_R2_001.fastq.gz

T1K_PATH=/opt/ohpc/Taiwania3/pkg/biology/T1K/T1K_v1.0.5

${T1K_PATH}/run-t1k \
-1 ${sampleR1} \
-2 ${sampleR2} \
--preset hla \
-f /work/${user}/T1K/hla_dna_seq.fa \
--od /work/${user}/HLA_typing \
-o S14
