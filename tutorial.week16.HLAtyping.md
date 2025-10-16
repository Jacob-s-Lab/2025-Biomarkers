# 20241218 Fantastic Genomic Biomarkers and Where to Find Them Practical Course (part X)
## The main content of this course

#### Learning to use `T1K` for HLA typing.

:warning:
### Prerequisite:Please first copy the necessary files for the course.
```markdown=
cd /work/username
mkdir T1K
cd T1K
rsync -avzP /work/u2499286/T1K/T1K_temp.sh ./
rsync -avzP /work/u2499286/hlaidx/hla_dna_seq.fa ./
```



ℹ️
#### What is T1K?
T1K (The ONE genotyper for KIR and HLA) is a computational tool designed to infer alleles of highly diverse genes, such as KIR and HLA. T1K leverages RNA-seq, WES, or WGS read alignment results to calculate allele abundance for the provided reference allele sequences and determines the actual alleles for each gene based on these abundance data.

T1K also offers subsequent analysis steps, including novel SNP (single nucleotide polymorphism) detection and single-cell expression analysis. T1K supports both single-end and paired-end sequencing data and is compatible with any read length.



### Step 1: Execute `T1K_temp.sh`
1. Open `T1K_temp.sh` and enter the correct path names.
    (1) <font color="#f00">Red</font> section: Replace the TA's username with your own.
    (2)<font color="#F7A004">Yellow</font>  section: Change it to your own file path.
    ![截圖 2024-12-17 下午1.40.53](https://hackmd.io/_uploads/H1bz190NJx.png)

       
2. Execute`T1K_temp.sh`
```
sbatch T1K_temp.sh
```
3. Result: A total of 7 files, all located in the HLA typing folder.
- S14_aligned_1.fa: Records the portion of the reference sequence that R1 aligned to.
- S14_aligned_2.fa: Records the portion of the reference sequence that R2 aligned to.
- S14_allele.tsv: Records the simplified HLA typing results.
- S14_allele.vcf: Records the variant sites identified by alignment to the HLA gene.
- S14_candidate_1.fq: R1 reads identified as likely related to HLA sequences.
- S14_candidate_2.fq: R2 reads identified as likely related to HLA sequences.
- S14_genotype.tsv: Records the main analysis results for reviewing HLA genotypes.
####  <font color="#f00">We need the HLA typing results for S14, which are all stored in the file `S14_genotype.tsv`.</font>
![截圖 2024-12-10 下午1.48.11](https://hackmd.io/_uploads/B1GNI8HEJe.png)
Column Descriptions (using the first row as an example):

- HLA-A: Gene name.
- 2: Number of gene alleles found.
- HLA-A*02:01:01: First gene allele.
- 90.834386: Confidence score.
- 35: Read depth.
- HLA-A*03:01:01: Second gene allele.
- 72.496481: Confidence score.
- 8: Read depth.
- HLA-A*03:279N: Additional possible gene allele.
- 16.559770: Score.
- 1: Read depth.






 ------------------
# 生物標記物與它們的產地實作課程(十)
## 本次課程主要內容
    
#### 學習利用 T1K 做 HLA typing。


:warning:
### 前情提要：請先複製課程所需檔案
```markdown=
cd /work/username
mkdir T1K
cd T1K
rsync -avzP /work/u2499286/T1K/T1K_temp.sh ./
rsync -avzP /work/u2499286/hlaidx/hla_dna_seq.fa ./
```


ℹ️
#### 什麼是 T1K？
T1K（The ONE genotyper for KIR and HLA）是一種計算工具，用於推測多樣性基因（如 KIR 和 HLA）的等位基因，T1K 基於 RNA-seq、WES 或 WGS 的讀取比對結果，針對提供的等位基因參考序列計算等位基因的豐富度，並使用這些豐富度數據來確定每個基因的真實等位基因。

T1K 還提供後續分析步驟，包括新穎 SNP（單核苷酸多樣性）的檢測和單細胞表現分析，T1K 支援單端（single-end）和雙端（paired-end）測序數據，且適用於任何讀取長度。


### step 1 : 執行 `T1K_temp.sh`

1. 打開 `T1K_temp.sh` 並寫入正確的路徑名稱。
    (1) <font color="#f00">紅色</font>部分：助教的 **username** 改成自己的。
    (2) <font color="#F7A004">黃色</font>部分：改成自己的檔案路徑。
![截圖 2024-12-17 下午1.40.53](https://hackmd.io/_uploads/HJm-1c0EJe.png)



    
2. 執行 `T1K_temp.sh`
```
sbatch T1K_temp.sh
```
3. 得到結果 : 共7個檔案,都在 HLA typing 資料夾中。
- S14_aligned_1.fa ：紀錄參考序列中R1比對到的部分。
- S14_aligned_2.fa ：紀錄參考序列中R2比對到的部分。
- S14_allele.tsv ：紀錄簡易的HLA分型結果。
- S14_allele.vcf ：紀錄與HLA基因比對的變異位點。
- S14_candidate_1.fq ：R1被認為與HLA相關的序列。
- S14_candidate_2.fq ：R2被認為與HLA相關的序列。
- S14_genotype.tsv ：紀錄主要分析結果用於查看HLA分型。

#### <font color="#f00">我們需要 S14 的 HLA typing 結果都在 `S14_genotype.tsv` 這個檔案中</font>。
![截圖 2024-12-10 下午1.48.11](https://hackmd.io/_uploads/B1GNI8HEJe.png)
欄位介紹（第一列為例）：

- HLA-A: 基因名稱。
- 2: 基因數量。
- HLA-A*02:01:01:第一個基因亞型。
- 90.834386：信賴分數。
- 35:深度。
- HLA-A*03:01:01：第二個基因亞型。
- 72.496481：信賴分數。
- 8:深度。
- HLA-A*03:279N：額外有可能的基因亞型。
- 16.559770：分數。
- 1:深度。
