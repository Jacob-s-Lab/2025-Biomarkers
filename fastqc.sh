#!/usr/bin/sh
#SBATCH -A ACD114093          # Account name/project number
#SBATCH -J fastqc                # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2                 # 使用的core數 請參考Queue資源設定
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.log   # Path to the standard output file
#SBATCH -e err.log
#SBATCH --mail-user=winche92221@gmail.com
#SBATCH --mail-type=END
# 國網使用

set -v -x
echo "start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"


# Please enter the R1 & R2 file name and your username
user=evelyn92
sampleR1=/work/${user}/result/fastqc/SRR13076392_S14_L002_R1_001.fastq.gz
sampleR2=/work/${user}/result/fastqc/SRR13076392_S14_L002_R2_001.fastq.gz
mkdir fastqc_S14

## Set up the environment for running fastqc
module load biology
module load FastQC


## Analyzing your sample's sequence QC by fastqc
echo "fastqc start"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

#fastqc -o ./fastqc_S14 ${sampleR1} ${sampleR2} 
fastqc ${sampleR1} ${sampleR2} -o fastqc_S14

echo "fastqc finished"
