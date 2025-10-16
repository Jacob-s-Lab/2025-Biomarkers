# 20241211 Fantastic Genomic Biomarkers and Where to Find Them Practical Course (part IX)
## The main content of this course


#### Learning to use `Manta` for SV calling.


:::warning
### Prerequisite: Please copy the required files for the course first. 
```markdown=
cd /work/username
mkdir SV
cd SV
rsync -avzP /work/u2499286/SV/format_change.py ./
rsync -avzP /work/u2499286/SV/HG002_SVs_Tier1_noVDJorXorY_v0.6.2_hs38DH.bed ./
rsync -avzP /work/u2499286/SV/manta_modified.sh ./
rsync -avzP /work/u2499286/SV/thalassemia_pipeline_course/demo/NGS1_20170103B.hs38DH.dedup.postalt.sorted.BQSR_chr16.bam ./
rsync -avzP /work/u2499286/SV/thalassemia_pipeline_course/demo/NGS1_20170103B.hs38DH.dedup.postalt.sorted.BQSR_chr16.bam.bai ./
```
:::

:::info
What is `Manta`?
- `Manta` is a tool for detecting genomic structural variations (SVs). The structural variants include large insertions, deletions, rearrangements, inversions, and other complex variation types. Illumina developed this tool, specifically designed for handling high-throughput sequencing (HTS) data.

What are **structural variants (SVs)**?
- **Structural variants** refer to DNA structural changes in chromosomes greater than 50 base pairs (bp). These variants may involve phenomena such as insertions, deletions, inversions, and duplications in chromosomes, potentially having profound impacts on gene function and the health of an organism.
:::

### Step 1: Execute `Manta`
1. Open `manta_modified.sh` and update the file with the correct path names.
    - Highlighted in <font color="#f00">red</font>: Replace the TA's username with your own username.
    ![截圖 2024-12-10 上午11.30.59](https://hackmd.io/_uploads/rkJuI4HV1l.png)
![截圖 2024-12-10 上午11.27.47](https://hackmd.io/_uploads/ByYd8ErEyg.png)

2. Execute `manta_modified.sh`
```
sbatch manta_modified.sh
```
3. Get two folders: `raw` and `filter`.
#### `raw` folder：
- candidateSV.vcf.gz
- candidateSV.vcf.gz.tbi
- candidateSmallIndels.vcf.gz
- candidateSmallIndels.vcf.gz.tbi
- diploidSV.vcf
- diploidSV.vcf.gz.tbi

#### `filter` folder：
#### <font color="#f00">The large deletion and insertion variants you want to check are all in this folder.</font>
**File names with "pass": Filtered for rows with the PASS field.**

**File names with "recode": Filtered for regions in the confident interval.**
- thalassemia_hg38_manta.vcf : All information.
- thalassemia_hg38_manta_del_pass.vcf
- thalassemia_hg38_manta_del_pass_filtered.recode.vcf
- thalassemia_hg38_manta_ins_pass.recode.vcf
- thalassemia_hg38_manta_ins_pass.vcf

### Step 2: Use IGV to confirm Deletion sites

<font color="#f00">This sample in the class is a thalassemia carrier. One large deletion can be observed in the positions of the *HBA1* and *HBA2* genes on chromosome 16 .</font>

:::warning
Command to open `IGV`:
```marksh
sh /opt/ohpc/Taiwania3/pkg/biology/IGV/IGV_v2.10.3/igv.sh
```
:::

1. Select input file:

    - Genome selection: hg38
    - File selection: NGS1_20170103B.hs38DH.dedup.postalt.sorted.BQSR_chr16.bam
    - Chr selection: chr16
    - Position selection: `chr16:162,130-189,352`
2. Right-click to select:
    - Expanded view
    - View as pairs
    - Color alignments by -> insert size and pair orientation
    - Sort alignments by -> insert size

![截圖 2024-12-05 下午3.34.43](https://hackmd.io/_uploads/B1kWGWx4Jg.png)




 ------------------------
# 生物標記物與它們的產地實作課程(九)
## 本次課程主要內容
    
#### 學習利用 `manta` 做 SV calling。


:::warning
### 前情提要：請先複製課程所需檔案
```markdown=
cd /work/username
mkdir SV
cd SV
rsync -avzP /work/u2499286/SV/format_change.py ./
rsync -avzP /work/u2499286/SV/HG002_SVs_Tier1_noVDJorXorY_v0.6.2_hs38DH.bed ./
rsync -avzP /work/u2499286/SV/manta_modified.sh ./
rsync -avzP /work/u2499286/SV/thalassemia_pipeline_course/demo/NGS1_20170103B.hs38DH.dedup.postalt.sorted.BQSR_chr16.bam ./
rsync -avzP /work/u2499286/SV/thalassemia_pipeline_course/demo/NGS1_20170103B.hs38DH.dedup.postalt.sorted.BQSR_chr16.bam.bai ./
```
:::

:::info
#### 什麼是 `manta` ？
- `manta`是一個用於檢測基因組結構變異（Structural Variations, SVs）的工具，處理的structural variants 包括大片段插入（Insertion）、缺失（Deletion）、重組（Rearrangement）、倒位（Inversion）以及其他複雜變異類型，此工具是由 Illumina 開發的，專為處理高通量測序（High-Throughput Sequencing, HTS）數據而設計。

#### 什麼是 structural variants (SVs)？
- 指相對於 reference genome ，染色體中大於50個 bp 的 DNA 結構改變。這些 variants 可能涉及染色體的插入、刪除、倒位、重複等現象，可能對基因功能和生物體健康產生深遠影響。


:::

### step 1 : 執行`manta`


1. 打開 `manta_modified.sh` 並寫入正確的路徑名稱。
     - <font color="#f00">紅色</font>部分：助教的 **username** 改成自己的。
    ![截圖 2024-12-10 上午11.30.59](https://hackmd.io/_uploads/HkLcIVBV1g.png)
![截圖 2024-12-10 上午11.27.47](https://hackmd.io/_uploads/r10qLVHV1e.png)



    
2. 執行 `manta_modified.sh`
```
sbatch manta_modified.sh
```
3. 得到兩個資料夾: `raw` 和 `filter`
#### `raw` 資料夾：
- candidateSV.vcf.gz
- candidateSV.vcf.gz.tbi
- candidateSmallIndels.vcf.gz
- candidateSmallIndels.vcf.gz.tbi
- diploidSV.vcf
- diploidSV.vcf.gz.tbi

#### `filter` 資料夾：
#### <font color="#f00">要看的大片段 deletion、insertion variants 都在這個資料夾</font>。
#### 檔名後有加上`pass`：篩選 filter 為 PASS 的欄位。
#### 檔名後有加上`recode`：篩選 confident interval 的區域。
- thalassemia_hg38_manta.vcf : 紀錄全部資訊。
- thalassemia_hg38_manta_del_pass.vcf
- thalassemia_hg38_manta_del_pass_filtered.recode.vcf
- thalassemia_hg38_manta_ins_pass.recode.vcf
- thalassemia_hg38_manta_ins_pass.vcf


### step 2 : 利用 `IGV` 確認 Deletion 位點
<font color="#f00">上課所使用的樣本是 thalassemia 的樣本，這種樣本會在 chr16 中 *HBA1* 及 *HBA2* genes 的位置看到大片段 deletion。</font>
:::warning
#### 打開 `IGV` 的指令：
```marksh
sh /opt/ohpc/Taiwania3/pkg/biology/IGV/IGV_v2.10.3/igv.sh
```

:::
1. 選擇輸入檔
    - Genome 選擇: hg38
    - File 選擇: NGS1_20170103B.hs38DH.dedup.postalt.sorted.BQSR_chr16.bam
    - chr 選擇: chr16
    - 位置選擇: `chr16:162,130-189,352`

2. 按右鍵選擇：
    - Expanded view
    - View as pairs
    - Color alignments by -> insert size and pair orientation
    - Sort alignments by -> insert size


![截圖 2024-12-05 下午3.34.43](https://hackmd.io/_uploads/B1kWGWx4Jg.png)
