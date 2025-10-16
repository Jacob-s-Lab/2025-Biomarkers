#!/usr/bin/sh
#SBATCH -A ACD114093        # Account name/project number
#SBATCH -J roc_HC          # Job name
#SBATCH -p ngscourse        # Partition name
#SBATCH -c 2               
#SBATCH --mem=13g           
#SBATCH -o roc_HC.out.log          # Path to the standard output file
#SBATCH -e roc_HC.err.log
#SBATCH --mail-user=
#SBATCH --mail-type=END

user=evelyn92
DIR_rop=/work/${user}/rocplot/rocplot_HC
mkdir -p ${DIR_rop}
cd ${DIR_rop}
Rscript /work/${user}/rocplot/rocplot_test.Rscript hap_plot -pr /work/${user}/hap/S18_HC_hap/output_prefix:S18 /work/${user}/hap/S15_HC_hap/output_prefix:S15 /work/${user}/hap/S14_HC_hap/output_prefix:S14
