#!/usr/bin/sh
#SBATCH -A ACD113120        # Account name/project number
#SBATCH -J fastqc         # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.log          # Path to the standard output file
#SBATCH -e err.log
#SBATCH --mail-user=
#SBATCH --mail-type=END

# Please enter the R1 & R2 file name and your username
R1=SRR13076392_S14_L002_R1_001
R2=SRR13076392_S14_L002_R2_001
user=u2499286
fastqdir=/work/u2499286/result/fastqc

R1_file=${fastqdir}/${R1}.fastq.gz
R2_file=${fastqdir}/${R2}.fastq.gz
## Set up the environment for running fastqc
module load biology
module load FastQC

## Analyzing your sample's sequence QC by fastqc
fastqc ${R1_file} ${R2_file} -o /work/${user}/result/
