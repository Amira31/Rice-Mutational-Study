#  **Manual to Rice Mutational Analysis**
This manual provides an end-to-end guide on running a mutational bioinformatics pipeline to discover candidate SNPs associated with the target mutant phenotype. This pipeline will be mainly executed with bash commands, Python and R scripts. However, users can freely choose which programming platforms they want to work in, at their convenience, depending on their work objectives. 

In this pipeline, I used Ubuntu 2404.68.0 (accessed with WSL) for data wrangling and manipulation, Python 3.12.3 for data visualisation, and Rstudio 2025.05.1 for data annotation and visualisation.

Before embarking on the analysis, users need to set up the appropriate environment and tools, if not already done so, for seamless execution. Users can read `here` on how to start setting up the environment. 


## Table of Contents
1. [Prerequisite Checklist](#step-0-prerequisite-checklist)
2. [Step 1: Setting up directory](#step-1-setting-up-directory)
3. [Step 2: Quality control and trimming](#step-2-quality-control-and-trimming)
4. [Step 3: Indexing of reference genome](#step-3-indexing-of-reference-genome)
5. [Step 4: Aligning FASTQ reads to FASTA reference](#step-4-aligning-fastq-reads-to-fasta-reference)
6. [Step 5: Duplicate marking and indexing of BAM](#step-5-duplicate-marking-and-indexing-of-bam)
7. [Step 6: Variant calling](#step-6-variant-calling)
8. [Step 7: Variant annotation](#step-7-variant-annotation)
   
## Step 0: Prerequisite Checklist

Before analysis, users must ensure that input data are provided as `FASTQ` files (either uncompressed or gzip-compressed) for both wild-type and mutant samples. 

`FASTQ` format is required because it has Phred quality scores for each base (in line 4) which is essential for reliable variant detection.

While assembled FASTA sequences can be used for sequence comparison, they lack base-quality and read-depth information and are therefore unsuitable for comprehensive variant analysis.



```bash
# Example of four-line string in FASTQ file
(line 1) @SEQ_ID
(line 2) GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTGA
(line 3) +
(line 4) !''*((((***+))%%%++)(%%%%).1***-+*''))**5555555555
```
Alternatively, users can also use test dataset available `here` in this Rice Mutational Study Github repository for testing. Once users have verified the pipeline's functionality and reproducibility, users may repeat with their real genomics dataset.  

## Step 1: Setting up directory

Users will need to create a directory to store their input and output datasets in one place. This is to ensure accessibility to these datasets when running bash commands. Bash commands can be tricky for beginners if files are in different places. 

**1. To create a directory:** `bash`

Users are suggested to create a parent directory and the subdirectories. Users may name the parent directory as the type of the analysis. For the sub-directories, users may name them according to the output format. 

```bash
# Get into the mounted D drive for WSL Linux:
cd /mnt
cd d

# Create parent & sub directories simultaneously
mkdir -p rice_wgrs/00_fastq
mkdir -p rice_wgrs/10_ref
mkdir -p rice_wgrs/20_bam
mkdir -p rice_wgrs/30_vcf
```
`mkdir` syntax creates a directory. `-p` syntax directs the failed attempt to be silent and return a normal status if there is already an existing directory.

These commands will give out folder structure as below:
```bash
rice_wgrs/
└── /00_fastq
└── /10_ref
└── /20_bam
└── /30_vcf 
```

* `rice_wgrs` Parent directory 
* `00_fastq` Subdirectory for storing FASTQ files
* `01_ref` Subdirectory for storing reference genome files and GFF
* `02_bam` Subdirectory for storing BAM files after alignment
* `03_vcf` Subdirectory for storing VCF files after variant calling

Additional subdirectories may be created to organise multiple datasets (e.g., FastQC reports for raw and trimmed reads), as long as their naming and placement are consistent with the directory structure shown above.

Example:
```bash
rice_wgrs/
└── /00_fastq
└── /01_trimmed_fastq
└── /02_fastqc
└── /03_trimmed_fastqc

└── /10_ref
└── /20_bam
└── /30_vcf 
```

## Step 2: Quality control and trimming

**1. Trimming reads** `bash`

As a general rule, users should process FASTQ reads prior to alignment analysis by removing adapter sequences, low-quality bases, and short reads to ensure data viability. The trimming and quality filtering can be performed using command-line tools such as `TRIMMOMATIC` or `fastp`. In this pipeline, I used  `TRIMMOMATIC v0.39`.

```bash
# create directory for trimmed FASTQ
mkdir -p 01_trimmed_fastq

# run trimming for wild-type:
trimmomatic PE -threads 8 \
  00_fastq/WT_R1.fastq.gz 00_fastq/WT_R2.fastq.gz \
  01_trimmed_fastq/WT_R1_paired.fq.gz 01_trimmed_fastq/WT_R1_unpaired.fq.gz \
  01_trimmed_fastq/WT_R2_paired.fq.gz 01_trimmed_fastq/WT_R2_unpaired.fq.gz \
  ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
  LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:75

# run trimming for mutant:
trimmomatic PE -threads 8 \
  00_fastq/M_R1.fastq.gz 00_fastq/M_R2.fastq.gz \
  01_trimmed_fastq/M_R1_paired.fq.gz 01_trimmed_fastq/M_R1_unpaired.fq.gz \
  01_trimmed_fastq/M_R2_paired.fq.gz 01_trimmed_fastq/M_R2_unpaired.fq.gz \
  ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
  LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:75

```

* `PE` Paired-end mode
* `-threads 8` Using 8 CPU threads
* `TruSeq3-PE-2.fa` Illumina adapter sequence type to remove
*  `LEADING:3` Removes bases with quality < 3 from the start of the read (5')
*  `TRAILING:3` Removes bases with quality < 3 from the end of the read (3')
*  `SLIDINGWINDOW:4:15` Sliding window of a size of 4 bases. The window slides from 5' to 3' and removes reads that have an average quality belows 15. 
*  `MINLEN:75` Removes reads that have < 75 bases

**2. Quality check** `bash`

After trimming low-quality bases and adapter sequences, the paired and unpaired FASTQ outputs should be assessed using `FastQC` to evaluate read quality. Here, I used `FastQC v0.12.1`. 

```bash
# create directory for FastQC report
mkdir -p 02_fastqc

# run FastQC analysis for wild-type
fastqc -t 8 \
  01_trimmed_fastq/WT_R1_paired.fq.gz \
  01_trimmed_fastq/WT_R2_paired.fq.gz \
  01_trimmed_fastq/WT_R1_unpaired.fq.gz \
  01_trimmed_fastq/WT_R2_unpaired.fq.gz \
  -o 02_fastqc

# run FastQC analysis for mutant
fastqc -t 8 \
  01_trimmed_fastq/M_R1_paired.fq.gz \
  01_trimmed_fastq/M_R2_paired.fq.gz \
  01_trimmed_fastq/M_R1_unpaired.fq.gz \
  01_trimmed_fastq/M_R2_unpaired.fq.gz \
  -o 02_fastqc
```
`fastqc` will assess the reads and provide a summary that consists of multiple modules:

* `Basic Statistics`
* `Per Base Sequence Quality`
* `Per Tile Sequence Quality`
* `Per Sequence Quality Scores`
* `Per Base Sequence Content`
* `Per Sequence GC Content`
* `Per Base N Content`
* `Sequence Length Distribution`
* `Sequence Duplication Levels`
* `Overrepresented Sequences`
* `Adapter Content`

Users are advised to run `FastQC` for both raw FASTQ and trimmed FASTQ, ensuring improved quality in each module for trimmed FASTQ, especially in terms of adapter content and per-base sequence quality. 

## Step 3: Indexing of reference genome 

First, users will need to download reference genome sequence in `FASTA` format `(.fa)`or `(.fna)` from NCBI together with its gene annotation file in GFF/GTF/GFF3 format. We will download these files with `wget` and unzip them with `gunzip`. 

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

Once users have their reference `.fna` or `.fa` file, they need to index the genome first. It is a rule of thumb to create `index` files for the reference genome before they can begin with the alignment process. 

Think of an index as a table of contents in a book. Searching for a specific subtopic from the table of contents page is faster than searching page by page from the beginning. That is similar to what reference indexing tries to achieve. An aligner tool will jump directly to specific read positions based on the index file, rather than browsing the genome from the beginning. 

Without an indexed reference, the aligner will work extremely slowly and often will fail. Even with index files, the aligner usually takes 2-10 hours to complete the entire ~370 Mb rice genome. 

To index the reference before alignment, either `BWA` Burrows-Wheeler Aligner, `BWA-MEM2` and `SAMtools` can be used for short-read WGS data. In this pipeline, I used `BWA v0.7.17`. 

In the future, users will need to view the aligned BAM on the Integrative Genomics Viewer (IGV) software. To be able to view the reference genome sequences, users are required beforehand to index the FASTA with `samtools faidx`. In this pipeline, I used `SAMtools v1.19.2`. 

```bash
# Index the reference with BWA for alignment  
bwa index 10_ref/Nipponbare.fna

#Index the reference with samtools for IGV viewing
samtools faidx 10_ref/Nipponbare.fna
```
Indexing will create these five index outputs:
* `Nipponbare.fna.amb`
* `Nipponbare.fna.ann`
* `Nipponbare.fna.bwt`
* `Nipponbare.fna.pac`
* `Nipponbare.fna.sa`
  
* `Nipponbare.fna.fai`

## Step 4: Aligning FASTQ reads to FASTA reference

For a mutational study, users need to align both wild-type and mutant read pairs (R1 and R2) separately to the reference genome. Direct alignment between non-reference wild-type and mutant is not a general standard of practice and can lead to false-positive variant calling. 

**1. Decide how many threads to use** `bash`

Before running the `bwa-mem` commands, users can check the CPU specifications of their laptop or PC to estimate how many threads they can allocate to the alignment. 

This is because alignment process is usually painstakingly long due to it being computationally intensive. So, activating multi-threading allows the aligner to split  different batches of reads into multiple parallel CPU queues to reduce the total runtime. Then, the aligner will merge back these parallel alignment records into a single SAM/BAM output stream.


```bash
# to check available CPU cores
nproc
```
* `16 logical CPU cores` Can use up to 16 threads.
* `8 logical CPU cores` Can use up to 8 threads.

**2. Alignment of short reads to reference genome** `bash`

```bash
# reads alignment to reference for wild-type
bwa mem -t 16 10_ref/Nipponbare.fna \
  01_trimmed_fastq/MR297_R1_paired.fq.gz \
  01_trimmed_fastq/MR297_R2_paired.fq.gz \
  > 20_bam/MR297.sam

# reads alignment to reference for mutant
bwa mem -t 16 10_ref/Nipponbare.fna \
  01_trimmed_fastq/ML-1_R1_paired.fq.gz \
  01_trimmed_fastq/ML-1_R2_paired.fq.gz \
  > 20_bam/ML-1.sam
```

* `mem` The algorithm used by BWA tool to align paired-end reads.
* `-t 16` Uses 16 CPU threads.
* `>` Directs the output to the specified SAM format.

**3. Conversion of SAM to unsorted BAM** `bash`

SAM is a huge sequence file (~40 GB). Using this file directly for the downstream process can be inefficient. Thus, converting SAM to unsorted BAM is highly recommended to reduce inefficiency. 

Please take note that at this step, BAM has not been sorted yet from query-based to coordinate-based. This is because during the duplicate marking step, `samtools fixmate` requires a query-based BAM input for it to fill MS (mate score) tags. 


```bash
# convert SAM to unsorted BAM
samtools view -@ 4 -b 20_bam/MR297.sam -o 20_bam/MR297_unsorted.bam
samtools view -@ 4 -b 20_bam/ML-1.sam -o 20_bam/ML-1_unsorted.bam
```
* `samtools view` Subcommand for manipulating SAM/BAM data
* `-@ 4` Uses 4 CPU threads
*  `-b` Converts the output into BAM format
* `-o` Specifies the output name
  

## Step 5: Duplicate marking and indexing of BAM

**1. Converting coordinate-sorted BAM to query-sorted BAM (optional):** `bash`

If you have performed `samtools sort` during the alignment process, it will produce an output of coordinate-sorted BAM (has SO:coordinate field). 

And if you still want to do PCR marking, you can use `samtools collate` to convert coordinate-sorted BAM into query name-sorted BAM (has SO:queryname field). 

This is because `samtools fixmate` cannot perform tag filling if it cannot find the QNAME (queryname) field on the BAM reads. 

```bash
# create a folder for temporary files
mkdir -p 20_bam/tmp/MR297_collate

# samtools collate on wild-type
samtools collate -@ 4 -o 20_bam/MR297_collated.bam 20_bam/MR297_sorted.bam -T 20_bam/tmp/MR297_collate

# samtools collate on mutant
samtools collate -@ 4 -o 20_bam/MR297_collated.bam 20_bam/MR297_sorted.bam -T 20_bam/tmp/MR297_collate
```
* `-@ 4` Uses 4 CPU threads
* `-o` Directs the collated output into new BAM
* ` -T 20_bam/tmp/MR297_collate` Direct the algorithm where to dump multiple temporary files 

**2. PCR duplicate marking** `bash`

In WGS sequencing, marking PCR duplicates is important to reduce false-positive reads. Oftentimes, sample preparation and PCR amplification before sequencing will lead to the same reads being replicated and sequenced. So, marking them with `samtools markdup` will tell variant caller tools to ignore these duplicated reads. 

However, to run `samtools markdup`, it is prerequisite to add mate tags to the reads with `samtools fixmate` and sort the reads by their coordinate position with `samtools sort` as below. 


```bash
# create directory for temporary files in Linux root
mkdir -p /root/tmp

# samtools fixmate on wild-type
samtools fixmate -m \
  20_bam/MR297_collated.bam \
  20_bam/MR297_fixmate.bam

# samtools fixmate on mutant
samtools fixmate -m \
  20_bam/ML-1_collated.bam \
  20_bam/ML-1_fixmate.bam

# samtools sort on wild-type
samtools sort -@ 2 -m 512M -T /root/tmp/MR297_sort \
  -o 20_bam/MR297_fixmate.sorted.bam \
  20_bam/MR297_fixmate.bam

# samtools sort on wild-type
samtools sort -@ 2 -m 512M -T /root/tmp/ML-1_sort \
  -o 20_bam/ML-1_fixmate.sorted.bam \
  20_bam/ML-1_fixmate.bam

# samtools markdup on wild-type
samtools markdup \
  20_bam/MR297_fixmate.sorted.bam \
  20_bam/MR297_markdup.bam

# samtools markdup on mutant
samtools markdup \
  20_bam/ML-1_fixmate.sorted.bam \
  20_bam/ML-1_markdup.bam
```
* `-m` Add and correct mate score tags and fields (MC, MQ, TLEN)
* `\` Tells the shell the command continues on the next line, thus, do not execute yet
* `@ 2` Uses 2 CPU threads
* `-m 512M` Working memory size per thread for sorting
* `-T` Directs samtools to dump temporary files under the specified file prefix
* `-o` Directs samtools where to write the stdout

These collective commands will result to multiple large intermediate BAM files (> 5GB each). So, it is advisable to remove the intermediate BAM files as below:
* `collated.bam`
* `fixmate.bam`
* `fixmate.sorted.bam`

```bash
# confirm file names
ls *_collated.bam *_fixmate.bam *_fixmate.sorted.bam

# remove unwanted intermediate BAM files
rm *_collated.bam && rm *_fixmate.bam && rm *_fixmate.sorted.bam
```

**3. Verifying BAM** `bash`
```bash
# verify the wild-type markdup BAM
samtools view -H 20_bam/MR297_markdup.bam | head
samtools view 20_bam/MR297_markdup.bam | head -n 5
samtools view -H 20_bam/MR297_markdup.bam | grep '^@HD'

# verify the mutant markdup BAM
samtools view -H 20_bam/ML-1_markdup.bam | head
samtools view 20_bam/ML-1_markdup.bam | head -n 5
samtools view -H 20_bam/ML-1_markdup.bam | grep '^@HD'
```

**4. Indexing BAM** `bash`
```bash
# index BAM
samtools index 20_bam/MR297_markdup.bam
samtools index 20_bam/ML-1_markdup.bam
```

* `|` BWA sends alignment output directly to samtools sort without creating intermediate SAM file
* `samtools sort -@ 8` coordinate-sorting in 8 parallel threads
* `-o` Sorted output will be written in BAM and stored in 20_bam folder

## Step 6: Variant calling

**1. Creating joint VCF** `bash`

```bash
# variant calling of wild-type and mutant
bcftools mpileup \
  -a AD,ADF,ADR,DP \
  -q 30 -Q 20 \
  -f 10_ref/Nipponbare.fna \
  20_bam/MR297_markdup.bam \
  20_bam/ML-1_markdup.bam | \
bcftools call -vm -f GQ,GP -O u | \
bcftools filter -i 'INFO/MQ>=40' -Oz -o 30_vcf/MR297_ML-1.vcf.gz
```

**2. Filter VCF for unique SNPs with SNP index ≥ 0.75** `bash`

```bash
# filter for unique SNPs with SNP index ≥ 0.75
bcftools view -v snps \
  -i 'FORMAT/AD[0:1]=0 && (FORMAT/AD[1:1]/(FORMAT/AD[1:0]+FORMAT/AD[1:1])>=0.75) && FORMAT/DP[0:0]>=10 && FORMAT/DP[1:0]>=10' \
  30_vcf/MR297_ML-1.vcf.gz \
  -Oz -o 30_vcf/ML-1_unique_snps_0.75.vcf.gz
```
**Filtering parameter:**
*`-v snps` to only extract SNP type variants 

Assuming, 
*`(FORMAT/AD)` = allele depth 
*`[0:1] / [1:1] / [1:0] / [0:0]` = [sample:allele]
*`[0]` = MR297
*`[1]` = ML-1
*`[0]` = REF
*`[1]` = ALT

Thus, 
*`FORMAT/AD[0:1]=0` to remove MR297 ALT SNPs 
*`(FORMAT/AD[1:1]/(FORMAT/AD[1:0]+FORMAT/AD[1:1]>=0.75)` SNP index = ALT reads / (REF reads + ALT reads) ≥ 0.75

`FORMAT/DP[1:0]>=10` to ensure MR297 has at least 10 reads at every position, otherwise absent also means no coverage


 **3. Indexing VCF for downstream uses** `bash`
 ```bash
tabix -p vcf 30_vcf/MR297_ML-1.vcf.gz
tabix -p vcf 30_vcf/ML-1_unique_snps_0.75.vcf.gz
```

## Step 7. Variant annotation
```bash
# create a directory for annotation output
mkdir -p 31_snpeff

# define path to the database
SNPEFF_HOME=/home/mira/miniconda3/envs/snpeff/share/snpeff-5.2-1

# verify the path
ls $SNPEFF_HOME/data/
ls $SNPEFF_HOME/data/genomes

# check FASTA header
grep "^>" $SNPEFF_HOME/data/genomes/Nipponbare.fa | head

# check GFF header
cut -f1 $SNPEFF_HOME/data/Nipponbare/genes.gff | sort | uniq

# check VCF header
bcftools view -H $SNPEFF_HOME/data/ML-1_unique_snps_0.75.vcf.gz | cut -f1 | sort -u

# convert from NC_089035.1 to Chr1 in VCF
nano $SNPEFF_HOME/data/vcf_rename.txt
bcftools annotate --rename-chrs $SNPEFF_HOME/data/vcf_rename.tx
t \
  -O z \
  -o $SNPEFF_HOME/data/ML-1_unique_snps_0.75_renamed.vcf.gz \
  $SNPEFF_HOME/data/ML-1_unique_snps_0.75.vcf.gz

# verify new VCF header
view -H $SNPEFF_HOME/data/ML-1_unique_snps_0.75_renamed.vcf.gz | cut -f1 | sort -u

# index new VCF
tabix -p vcf $SNPEFF_HOME/data/ML-1_unique_snps_0.75_renamed.vcf.gz
ls $SNPEFF_HOME/data/ML-1_unique_snps_0.75_renamed.vcf.gz.tbi

# run snpeff
java -Xmx6g -jar $SNPEFF_HOME/snpEff.jar \
  -v Nipponbare \
  -stats 33_snpeff/ML-1_snpeff_summary.html \
  $SNPEFF_HOME/data/ML-1_unique_snps_0.75_renamed.vcf.gz \
  > 33_snpeff/ML-1_unique_snps_0.75_renamed.ann.vcf.gz

# convert annotated VCF to tab-delimited TXT
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' \
  33_snpeff/ML-1_unique_snps_0.75_renamed.ann.vcf.gz | \
awk -F'\t' 'BEGIN{OFS="\t"}{
    split($5,a,",");
    split(a[1],b,"|");
    print $1,$2,$3,$4, \
          b[2],b[3],b[4],b[5],b[6],b[7],b[8],b[9]

}' > 33_snpeff/ML-1_snpeff_ann_table.txt
```




## Mapping Statistics
### Substep: To calculate SNs substitution counts
```bash
bcftools query -l 30_vcf/MR297_ML-1.vcf.gz

bcftools view -s '20_bam/MR297_markdup.bam' -v snps -m2 -M2 30_vcf/MR297_ML-1.vcf.gz \
| bcftools view -i 'GT="0/1" || GT="1/1"' \
| bcftools query -f '%REF\t%ALT\n' \
| awk '{print $1 ">" $2}' \
| sort | uniq -c | sort -k2

bcftools view -s '20_bam/ML-1_markdup.bam' -v snps -m2 -M2 30
_vcf/MR297_ML-1.vcf.gz \
| bcftools view -i 'GT="0/1" || GT="1/1"' \
| bcftools query -f '%REF\t%ALT\n' \
| awk '{print $1 ">" $2}' \
| sort | uniq -c | sort -k2
```





##
