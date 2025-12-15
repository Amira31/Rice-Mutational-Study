#  **Manual to Rice Mutational Analysis**
This manual provides an end-to-end guide on running a mutational bioinformatics pipeline to discover candidate SNPs associated with the target mutant phenotype. This pipeline will be mainly executed with bash commands, Python and R scripts. However, users can freely choose which programming platforms they want to work in, at their convenience, depending on their work objectives. 

In this pipeline, I used Ubuntu 2404.68.0 (accessed with WSL) for data wrangling and manipulation, Python 3.12.3 for data visualisation, and Rstudio 2025.05.1 for data annotation and visualisation.

Before embarking on the analysis, users need to set up the appropriate environment and tools, if not already done so, for seamless execution. Users can read `here` on how to start setting up the environment. 


## Table of Contents
1. [Prerequisite Checklist](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
2. [Step 1: Setting up the directory](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
3. [Step 2: Quality control and trimming](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
4. [Step 3: Indexing of reference genome](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
5. [Prerequisite Checklist](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
6. [Setting up the directory](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
7. [Indexing of reference genome](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
8. [Prerequisite Checklist](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
9. [Setting up the directory](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
   
## Step 0: Prerequisite Checklist
Before starting, users need to ensure that the starting input file is in the `FASTQ` or zipped `FASTQ` format for their wild-type and mutant sequences. 

Users can also, optionally, use raw `FASTA` wild-type/mutant sequences. However, FASTA will not provide enough data for a comprehensive variant analysis. FASTQ, on the other hand, has a Phred quality score for each base (in line 4) thus allowing more depth to the analysis.

```bash
# Example of four-line string in FASTQ file
(line 1) @SEQ_ID
(line 2) GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTGA
(line 3) +
(line 4) !''*((((***+))%%%++)(%%%%).1***-+*''))**5555555555
```
Alternatively, users can also use test dataset available `here` in this Rice Mutational Study Github repository for testing. Once users have gone through the bottom and verified the pipeline's functionality, users may repeat with their real genomics dataset.  

## Step 1: Setting up directory

Users will need to create a directory to store their input and output datasets in one place. This is to ensure accessibility to these datasets when running bash commands. Bash commands can be tricky for beginners if files are in different places. 

**1. To create a directory:** `bash`

I would suggest users to create a parent directory and the sub-directories. Users may name the parent directory as the type of the analysis. For the sub-directories, users may name them according to the output format. 

```bash
# Get into the mounted D drive:
cd /mnt
cd d

# Create parent & sub- directories simultaneously
mkdir -p rice_wgrs/00_fastq
mkdir -p rice_wgrs/10_ref
mkdir -p rice_wgrs/20_bam
mkdir -p rice_wgrs/30_vcf
```

This will give out folder structure as below:
```bash
rice_wgrs/
└── /00_fastq
└── /10_ref
└── /20_bam
└── /30_vcf 
```
* `rice_wgrs` Parent directory 
* `00_fastq` For storing filtered FASTQ files
* `01_ref` For storing reference genome files and GFF
* `02_bam` For storing BAM files after alignment
* `03_vcf` For storing VCF files after variant calling

## Step 2: Quality control and trimming

**1. Trimming reads** `bash`

Before using FASTQ reads for alignment, users should process them to remove adapter sequences, low-quality bases, and reads that are too short. This trimming and quality filtering can be performed using tools such as `TRIMMOMATIC` or `fastp`. In this pipeline, I used `TRIMMOMATIC v0.39`. 

```bash
# create directory for trimmmed FASTQ
mkdir -p 01_trimmed_fastq

# run trimming for wild-type:
trimmomatic PE -threads 8 \
  00_fastq/MR297_R1.fastq.gz 00_fastq/MR297_R2.fastq.gz \
  01_trimmed_fastq/MR297_R1_paired.fq.gz 01_trimmed_fastq/MR297_R1_unpaired.fq.gz \
  01_trimmed_fastq/MR297_R2_paired.fq.gz 01_trimmed_fastq/MR297_R2_unpaired.fq.gz \
  ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
  LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:75

# run trimming for mutant:
trimmomatic PE -threads 8 \
  00_fastq/ML-1_R1.fastq.gz 00_fastq/ML-1_R2.fastq.gz \
  01_trimmed_fastq/ML-1_R1_paired.fq.gz 01_trimmed_fastq/ML-1_R1_unpaired.fq.gz \
  01_trimmed_fastq/ML-1_R2_paired.fq.gz 01_trimmed_fastq/ML-1_R2_unpaired.fq.gz \
  ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
  LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:75

```

* `PE` Paired-end mode
* `-threads 8` Using 8 CPU threads
* `TruSeq3-PE-2.fa` Illumina adapter sequences to remove
*  `LEADING:3` Removes bases with quality < 3 from the start of the read (5')
*  `TRAILING:3` Removes bases with quality < 3 from the end of the read (3')
*  `SLIDINGWINDOW:4:15` Sliding window of a size of 4 bases. The window slides from 5' to 3' and removes reads that have an average quality belows 15. 
*  `MINLEN:75` Removes reads that are < 75 bases

**2. Quality check** `bash`

After trimming low-quality bases and adapter sequences, the paired and unpaired FASTQ outputs should be assessed using `FastQC` to evaluate read quality. Here, I used `FastQC v0.12.1`. 

```bash
# create directory for FastQC report
mkdir -p 02_fastqc

# run FastQC analysis for wild-type
fastqc -t 8 \
  01_trimmed_fastq/MR297_R1_paired.fq.gz \
  01_trimmed_fastq/MR297_R2_paired.fq.gz \
  01_trimmed_fastq/MR297_R1_unpaired.fq.gz \
  01_trimmed_fastq/MR297_R2_unpaired.fq.gz \
  -o 02_fastqc

# run FastQC analysis for mutant
fastqc -t 8 \
  01_trimmed_fastq/ML-1_R1_paired.fq.gz \
  01_trimmed_fastq/ML-1_R2_paired.fq.gz \
  01_trimmed_fastq/ML-1_R1_unpaired.fq.gz \
  01_trimmed_fastq/ML-1_R2_unpaired.fq.gz \
  -o 02_fastqc
```

## Step 3: Indexing of reference genome 

First, users will need to download reference genome sequence in `FASTA` format `(.fa)` `(.fna)` from NCBI together with its gene annotation file in GFF/GTF/GFF3 format. We will download these files with `wget` and unzip them with `gunzip`. 

For this pipeline, I will use a japonica cv. Nipponbare AGIS 1.0 reference genome. 

**1. Retrieve FASTA and GFF links:** `bash`

Users might wonder how one can obtain independent links for FASTA and GFF because NCBI genome assembly main page only shows the bulk download menu. Fret not, all they need to do is click the `FTP` menu and it will redirect to a page with a list of hyperlinks. Browse through and look for links that end with `.fna.gz` and `gff.gz`. Now, users need to copy the `ftp` page link and add behind it the FASTA and GFF links they found on the terminal, as below:

```bash
# Download reference FASTA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.fna.gz   

# Download reference GFF
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.gff.gz
```
Since the file names are full of numbers, I would prefer to rename them into readable names.

```bash
# Rename FASTA
mv GCF_034140825.1_ASM3414082v1_genomic.fna.gz 10_ref/Nipponbare.fna.gz

# Rename GFF
mv GCF_034140825.1_ASM3414082v1_genomic.gff.gz 10_ref/Nipponbare.gff.gz
```
Now, users need to decompress these files so they become usable for shell scripting.
```bash
gunzip 10_ref/Nipponbare.fna.gz
gunzip 10_ref/Nipponbare.gff.gz
```

**2. Index the reference FASTA** `bash`

Now, once users have their reference `.fna` `.fasta` file, they need to index the genome first. It is a rule of thumb to create `index` files for the reference genome before they can begin with the alignment process. 

Think of an index as a table of contents in a book. Searching for a specific subtopic from the table of contents page is faster than searching page by page from the beginning. That is similar to what reference indexing tries to achieve. An aligner tool will jump directly to specific read positions based on the index file, rather than browsing the genome from the beginning. 

Without an indexed reference, the aligner will work extremely slowly and often will fail. Even with index files, the aligner usually takes 2-8 hours to complete the entire ~370 Mb rice genome. 

To index the reference, either `BWA` Burrows-Wheeler Aligner, `BWA-MEM2` and `SAMtools` can be used for short-read WGS data. In this pipeline, I used `BWA v0.7.17`. 

```bash
# Index the reference using BWA  
bwa index 10_ref/Nipponbare.fna
```
Indexing will create these five index outputs:
* `Nipponbare.fna.amb`
* `Nipponbare.fna.ann`
* `Nipponbare.fna.bwt`
* `Nipponbare.fna.pac`
* `Nipponbare.fna.sa`

## Step 4: Read alignment

Since they are doing a mutational study, users need to align both wild-type and mutant read pairs (R1 and R2) separately to the reference genome. Direct alignment between non-reference wild-type and mutant is not a general standard of practice and can lead to false-positive variant calling. 

**1. Decide how many threads to use** `bash`

Before running the `bwa-mem` or `bwa-mem2 mem` commands, users can check the CPU specifications of their laptop or PC to estimate how many threads they can allocate to the alignment. 

This is because alignment process is usually painstakingly long due to it being computationally intensive. So, activating multi-threading allows the aligner to split  different batches of reads into parallel queues of multiple CPU threads to reduce the total runtime. Then, the aligner will merge back these parallel alignment records into a single SAM/BAM output stream.


```bash
# to check available CPU cores
nproc
```
* `16 logical CPU cores` recommended to use 2-16 threads
* `8 logical CPU cores` recommended to use 2-8 threads

**2. Alignment of short reads to reference genome** `bash`

```bash
# reads alignment to reference for wild-type
bwa mem -t 4 10_ref/Nipponbare.fna \
  01_trimmed_fastq/MR297_R1_paired.fq.gz \
  01_trimmed_fastq/MR297_R2_paired.fq.gz | \
samtools sort -@ 8 -o 20_bam/MR297.sorted.bam

# reads alignment to reference for mutant
bwa mem -t 4 10_ref/Nipponbare.fna \
  01_trimmed_fastq/ML-1_R1_paired.fq.gz \
  01_trimmed_fastq/ML-1_R2_paired.fq.gz | \
samtools sort -@ 8 -o 20_bam/ML-1.sorted.bam
```
```bash
# 1) Alignment only (fast, low memory pressure)
bwa mem -t 4 10_ref/Nipponbare.fna \
  01_trimmed_fastq/MR297_R1_paired.fq.gz \
  01_trimmed_fastq/MR297_R2_paired.fq.gz \
  > 20_bam/MR297.sam

# 2) Sort + compress afterwards (controlled memory)
samtools sort \
  -@ 2 \
  -m 512M \
  -T ~/wgrs_tmp/MR297 \
  -o 20_bam/MR297.sorted.bam \
  20_bam/MR297.sam

# 3) Cleanup
rm 20_bam/MR297.sam
```



* `bwa mem -t 4` aligner tool that processes reads in 4 parallel threads
* `|` BWA sends alignment output directly to samtools sort without creating intermediate SAM file
* `samtools sort -@ 8` coordinate-sorting in 8 parallel threads
* `-o` Sorted output will be written in BAM and stored in 20_bam folder

## Step 5: BAM duplicate marking and indexing

```bash
# duplicate marking and indexing of wild-type BAM
samtools markdup -r \
  20_bam/MR297.sorted.bam \
  20_bam/MR297.bam

samtools index 20_bam/MR297.bam

# duplicate marking and indexing of mutant BAM
samtools markdup -r \
  20_bam/ML-1.sorted.bam \
  20_bam/ML-1.bam

samtools index 20_bam/ML-1.bam
```
* `samtools markdup`

## Step 6: Variant calling

```bash
# variant calling of wild-type
bcftools mpileup -a AD,ADF,ADR -B -q 30 -Q 20 -C 50 -f Nipponbare.fna \
20_bam/MR297.bam | \
bcftools call -vm -f GQ,GP -O u | \
bcftools filter -i 'INFO/MQ>=40 && INFO/DP>=10 && INFO/DP<=200' -O z -o 30_vcf/MR297.vcf.gz

# variant calling of mutant
bcftools mpileup -a AD,ADF,ADR -B -q 30 -Q 20 -C 50 -f Nipponbare.fna \
20_bam/ML-1.bam | \
bcftools call -vm -f GQ,GP -O u | \
bcftools filter -i 'INFO/MQ>=40 && INFO/DP>=10 && INFO/DP<=200' -O z -o 30_vcf/ML1.vcf.gz
```
Or, alternatively,

```bash
# join VCF
bcftools mpileup -a AD,ADF,ADR -B -q 30 -Q 20 -C 50 -f Nipponbare.fna \
20_bam/MR297.bam 20_bam/ML-1.bam | \
bcftools call -vm -f GQ,GP -O u | \
bcftools filter -i 'INFO/MQ>=40 && INFO/DP>=10 && INFO/DP<=200' -O z -o 30_vcf/MR297_ML-1.vcf.gz
```

```bash
tabix -p vcf 30_vcf/MR297_ML1.vcf.gz
```

```bash
# extract unique SNPs
bcftools view -v snps -i 'GT[ML-1]!=0 && GT[MR297]==0' \
30_vcf/MR297_ML1.vcf.gz -Oz -o 30_vcf/ML-1_unique_snps.vcf.gz

# add AF column with +fill-tags
bcftools +fill-tags ML-1_unique_snps.vcf.gz -- -t AF -Oz -o ML-1_unique_snps_AF.vcf.gz
```

optional:
```
# verify that AF column is there
bcftools view -h ML-1_unique_snps_AF.vcf.gz | grep "^##FORMAT=<ID=AF"

# see the values in the AF column
bcftools query -f '%CHROM\t%POS\t[%AF]\n' ML1_unique_snps_AF.vcf.gz | head
```
```
# filter for AF > 0.75
bcftools view -i 'AF[ML-1]>0.75' ML1_unique_snps_AF.vcf.gz -Oz -o ML-1_unique_snps_0.75.vcf.gz
```








<br/>
<br/>
<br/>

--- 

If users use the classic `BWA`:

```bash
# Align the wild-type to reference genome
bwa mem -t 16 rice_wgrs/10_ref/Nippombare.fna rice_wgrs/WT_R1.fq.gz WT_R2.fq.gz > WT.sam

# Align the mutant to reference genome
bwa mem -t 16 Nippombare.fna M_R1.fq.gz M_R2.fq.gz > M.sam
```
If users use `BWA-MEM2`:
```bash
# Align the wild-type to reference genome
bwa-mem2 mem -t 16 Nippombare.fna WT_R1.fq.gz WT_R2.fq.gz > WT.sam

#  Align the mutant to reference genome
bwa-mem2 mem -t 16 Nippombare.fna M_R1.fq.gz M_R2.fq.gz > M.sam
```

```bash
bwa mem -t 16 10_ref/Nipponbare.fna \
  01_trimmed_fastq/MR297_R1_paired.fq.gz \
  01_trimmed_fastq/MR297_R2_paired.fq.gz | \
samtools sort -@ 8 -o 20_bam/MR297.sorted.bam

samtools markdup -r \
  20_bam/MR297.sorted.bam \
  20_bam/MR297.markdup.bam

samtools index 20_bam/MR297.markdup.bam
```













## Substep: Filter VCF for SNP index = 1
Next, you will need to screen the large number of unique SNPs in the VCF to retain only SNPs with SNP index equals 1 (=1), or more than 0.8 (>= 0.8) if you want more selection. Then, optionally, you can validate the SNPs you filtered by visualising the BAM mutant/wild-type sequences on the IGV software. 

**1. Create a zip file for the VCF** `bash` 
 ```bash
bgzip MR297.vcf
bgzip ML-1.vcf
```

 This will give the output to MR297.vcf.gz and ML-1.vcf.gz

</br> **2. Check the VCF header**
 </br> Before filtering, you need to make sure SNP_INDEX ==column== name is indeed in the header list of the VCF. 

</br> **3. Set file names** `bash`

wild-type:
```bash
VCF_IN = "MR297.vcf.gz"
ANNOT = "MR297.annot"
HEADER = "MR297.header"
VCF_OUT = "MR297_allindex.vcf.gz"
VCF_INDEX1 = "MR297_index1.vcf.gz"
```

mutant:
```
VCF_INm = "ML-1.vcf.gz"
ANNOTm = "ML-1.annot"
HEADERm = "ML-1.header"
VCF_OUTm = "ML-1_allindex.vcf.gz"
VCF_INDEX1m = "ML-1_index1.vcf.gz"
```
<br/> **3. Create ANNOT file** `bash`

wild-type:
```
bcftools query -f '%CHROM\t%POS\t%DP4\n' $VCF_IN \
  | awk 'BEGIN {OFS = "\t"}{
      split($3, a, ",");
      ref = a[1] + a[2];
      alt = a[3] + a[4];
      if (ref + alt > 0) {
          snp_index = alt / (ref + alt);
      } else {
          snp_index = ".";
      }
      print $1, $2, snp_index;
  }' > "$ANNOT"
```

mutant:
```
bcftools query -f '%CHROM\t%POS\t%DP4\n' $VCF_INm \
  | awk 'BEGIN {OFS = "\t"}{
      split($3, a, ",");
      ref = a[1] + a[2];
      alt = a[3] + a[4];
      if (ref + alt > 0) {
          snp_index = alt / (ref + alt);
      } else {
          snp_index = ".";
      }
      print $1, $2, snp_index;
  }' > "$ANNOTm"
```

<br/> **4. Compress and index ANNOT** `bash`

wild-type:
```
bgzip -c $ANNOT > ${ANNOT}.gz
tabix -s1 -b2 -e2 ${ANNOT}.gz
```

mutant:
```
bgzip -c $ANNOTm > ${ANNOTm}.gz
tabix -s1 -b2 -e2 ${ANNOTm}.gz
```

<br/> **5. Create header for SNP_INDEX** `bash`

wild-type:
```
echo '##INFO=<ID=SNP_INDEX,Number=1,Type=Float,Description="SNP index calculated from DP4">' > $HEADER
```

mutant:
```
echo '##INFO=<ID=SNP_INDEX,Number=1,Type=Float,Description="SNP index calculated from DP4">' > $HEADERm
```

<br/> **6. Annotate input VCF with SNP_INDEX** `bash`

wild-type:
```
bcftools annotate \
  -a ${ANNOT}.gz \
  -h $HEADER \
  -c CHROM, POS, SNP_INDEX \
  -Oz -o $VCF_OUT \
  $VCF_IN

bcftools index $VCF_OUT
```

mutant:
```
bcftools annotate \
  -a ${ANNOTm}.gz \
  -h $HEADER \
  -c CHROM, POS, SNP_INDEX \
  -Oz -o $VCF_OUTm \
  $VCF_INm

bcftools index $VCF_OUTm
```
<br/>**7. Filter VCF_OUT for SNP Index = 1** `bash`

wild-type:
```
bcftools filter \
  -i "INFO/SNP_INDEX == 1" \
  $VCF_OUT \
  -Oz -o $VCF_INDEX1

bcftools index $VCF_INDEX1
```

mutant:
```
bcftools filter \
  -i "INFO/SNP_INDEX == 1" \
  $VCF_OUTm \
  -Oz -o $VCF_INDEX1m

bcftools index $VCF_INDEX1m
```
<br/>**8. Verify:** `bash`

wild-type:
```
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%SNP_INDEX\n' $VCF_INDEX1 | head
```

mutant:
```
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%SNP_INDEX\n' $VCF_INDEX1m | head
```

expected result printed on terminal:
```
Chr1    35477   A       G       1
Chr1    35884   C       T       1
Chr1    38073   A       G       1
Chr1    38463   C       T       1
Chr1    42795   T       C       1
Chr1    46296   T       C       1
Chr1    50877   A       G       1
Chr1    71547   A       G       1
Chr1    85087   T       G       1
Chr1    85349   C       T       1
```

```bash
# Check mapping stats
bcftools stats  MR219.vcf.gz
```

```bash
# CHeck contig names in GFF
cut -f1 01_ref/Nipponbare.gff | sort | uniq
```

```bash
# activate conda environment
conda activate bioinfo
```

```bash
gunzip -c 10_ref/Nipponbare1.fna.gz \
  | grep "^>" \
  | sed 's/^>//' \
  | cut -d' ' -f1 \
  | sort -u
```








##
