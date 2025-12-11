#  **Manual to Rice Mutational Analysis**
This manual provides an end-to-end guide on running a mutational bioinformatics pipeline to discover candidate SNPs associated with the target mutant phenotype. This pipeline will be mainly executed with bash commands, Python and R scripts. However, users can freely choose which programming platforms they want to work in, at their convenience, depending on their work objectives. 

In this pipeline, I used Ubuntu 2404.68.0 (accessed with WSL) for data wrangling and manipulation, Python 3.12.3 for data visualisation, and Rstudio 2025.05.1 for data annotation and visualisation.

Before embarking on the analysis, users need to set up the appropriate environment and tools, if not already done so, for seamless execution. Users can read `here` on how to start setting up the environment. 


## Table of Contents
1. [Prerequisite Checklist](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
2. [Setting up the directory](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
3. [Indexing of reference genome](https://github.com/Amira31/Rice-Mutational-Study/edit/main/01_Manual.md)
4. d
5. e
6. f
7. g
8. h
   
## Step 0: Prerequisite Checklist
Before starting, users need to ensure that the starting input file is in the `FASTQ` `(.fastq/.fq)` or zipped `FASTQ` `(.fastq.gz/.fq.gz)` format for their wild-type and mutant sequences. 

Users can also, optionally, use raw `FASTA` `(.fasta/.fa)` wild-type/mutant sequences. However, FASTA will not provide enough data for a comprehensive variant analysis. FASTQ, on the other hand, has a Phred quality score for each base (in line 4) thus allowing more depth to the analysis.

```
# Example of four-line string in FASTQ file.

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

```
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
```
rice_wgrs/
└── /00_fastq
└── /10_ref
└── /20_bam
└── /30_vcf 
```
* **rice_wgrs:** fff
* **00_fastq:**
* **01_ref:**
* **02_bam:**
* **03_vcf:**

## Step 2: Indexing of reference genome 

First, you will need to download reference genome sequence in `FASTA` format `(.fa/.fna)` from NCBI together with its gene annotation file in GFF/GTF/GFF3 format. We will download these files with `wget` and unzip them with `gunzip`. 

For this pipeline, I will use a japonica cv. Nipponbare AGIS 1.0 reference genome. 

**1. Retrieve FASTA and GFF links:** `bash`

You might wonder how one can obtain independent links for FASTA and GFF because NCBI genome assembly main page only shows the bulk download menu. Fret not, all you need to do is click the `FTP` menu and it will redirect to a page with a list of hyperlinks. Browse through and look for links that end with `.fna.gz` and `gff.gz`. Now, you need to copy the `ftp` page link and add behind it the FASTA and GFF links you found on the terminal, as below:

```
# Download reference FASTA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.fna.gz   

# Download reference GFF
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.gff.gz
```
Since the file names are full of numbers, I would prefer to rename them into readable names.

```
# Rename FASTA
mv GCF_034140825.1_ASM3414082v1_genomic.fna.gz /10_ref/Nipponbare.fna.gz

# Rename GFF
mv GCF_034140825.1_ASM3414082v1_genomic.gff.gz /10_ref/Nipponbare.gff.gz
```
Now, you need to decompress these files so they become usable for shell scripting.
```
gunzip /10_ref/Nipponbare.fna.gz
gunzip /10_ref/Nipponbare.gff.gz
```

**2. Index the reference FASTA** `bash`

Now, once you have your reference `.fna/.fasta` file, you need to index the genome first. It is a rule of thumb to create an `index` file for the reference genome before you begin with the alignment process. 

Think of an index as a table of contents in a book. Searching for a specific subtopic from the table of contents page is faster than searching page by page from the beginning. That is similar to what reference indexing tries to achieve. An aligner tool will jump directly to specific read positions based on the index file, rather than browsing the genome from the beginning. 

Without an index file, your aligner will work extremely slowly and often will fail. Even with an index file, an aligner usually takes 2-8 hours to complete the entire ~370 Mb rice genome. 

To index the reference, either `BWA` Burrows-Wheeler Aligner and `SAMtools` can be used for short-read WGS data. Certainly, you will need to install these tools first with `apt`. 

```
# Index the reference using BWA  
bwa index /10_ref/Nipponbare.fna

#Index the reference using SAMtools
samtools faidx /10_ref/Nipponbare.fna
```


## Step 3: Alignment of short reads to reference genome

Since you're doing a mutational study, you need to align both wild-type and mutant read pairs (R1 and R2) separately to the reference genome.  


“Before running the `bwa-mem` command, check your CPU specifications so you can estimate how many threads you can allocate for the alignment. 

If you're using a HPC node


```
bwa mem -t 16 /10_ref/Nipponbare.fna 
```

## Substep: Filter VCF for SNP index = 1
Next, you will need to screen the large number of unique SNPs in the VCF to retain only SNPs with SNP index equals 1 (=1), or more than 0.8 (>= 0.8) if you want more selection. Then, optionally, you can validate the SNPs you filtered by visualising the BAM mutant/wild-type sequences on the IGV software. 

**1. Create a zip file for the VCF** `bash` 
 ```
bgzip MR297.vcf
bgzip ML-1.vcf
```

 This will give the output to MR297.vcf.gz and ML-1.vcf.gz

</br> **2. Check the VCF header**
 </br> Before filtering, you need to make sure SNP_INDEX ==column== name is indeed in the header list of the VCF. 

</br> **3. Set file names** `bash`

wild-type:
```
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

















##
