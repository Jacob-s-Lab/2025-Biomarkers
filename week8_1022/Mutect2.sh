gatk Mutect2 \
  -R /opt/ohpc/Taiwania3/pkg/biology/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta \
  -I SRR13076392_S14_L002_.sorted.markdup.bam \
  -O SRR13076392_S14_L002_.M2.vcf.gz \
  -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9  \
  -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17  \
  -L chr18 -L chr19 -L chr20 -L chr21 -L chr22


echo "$(date '+%Y-%m-%d %H:%M:%S') Job finished"
