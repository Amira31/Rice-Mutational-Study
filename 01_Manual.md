#  **Manual to Rice Mutational Analysis**
This manual provides an end-to-end guide on running a mutational bioinformatics pipeline to discover candidate SNPs associated with the target mutant phenotype. This pipeline will be mainly executed with bash commands, Python and R scripts. However, users can freely choose which platforms they want to work in, at their convenience, depending on their work objectives. 

In this pipeline, I used Ubuntu 2404.68.0 (accessed with WSL) for data wrangling and manipulation, Python 3.12.3 for data visualisation, and Rstudio 2025.05.1 for data annotation and visualisation.

Before embarking on the analysis, users need to set up the appropriate environment and tools, if not already done so, for seamless execution. Users can read more `here` on how to start setting up the environment. 


## Table of Contents

## Step 0: Prerequisite Checklist


## Step 1: Alignment of short reads to reference genome

```
bwa index Nipponbare.fna
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

















##
