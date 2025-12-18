
# Get into the mounted D drive for WSL Linux:
cd /mnt
cd d

# Create parent & sub directories simultaneously
mkdir -p rice_wgrs/00_fastq
mkdir -p rice_wgrs/10_ref
mkdir -p rice_wgrs/20_bam
mkdir -p rice_wgrs/30_vcf

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


# Download reference FASTA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.fna.gz   

# Download reference GFF
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.gff.gz

# Rename reference FASTA
mv GCF_034140825.1_ASM3414082v1_genomic.fna.gz 10_ref/Nipponbare.fna.gz

# Rename reference GFF
mv GCF_034140825.1_ASM3414082v1_genomic.gff.gz 10_ref/Nipponbare.gff.gz

# decompress reference FASTA
gunzip 10_ref/Nipponbare.fna.gz
gunzip 10_ref/Nipponbare.gff.gz













