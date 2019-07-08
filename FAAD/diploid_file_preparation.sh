#!/bin/bash -e


source bcftools-1.2

srun bcftools view   -S dinaric.txt  samples.ann.vcf.gz  > dinaric.vcf
srun bcftools view   -S baltic.txt  samples.ann.vcf.gz  > baltic.vcf
srun bcftools view   -S se_carpathian.txt  samples.ann.vcf.gz  > se_carpathian.vcf
srun bcftools view   -S pannonean.txt  samples.ann.vcf.gz  > pannonean.vcf
srun bcftools view   -S w_carpathian.txt  samples.ann.vcf.gz  > w_carpathian.vcf

source snpEFF-4.3 
source jre-7.21 




srun java -jar  /nbi/software/testing/snpEFF/4.3/x86_64/SnpSift.jar extractFields dinaric.missense.vcf CHROM POS REF ALT AC AN  "ANN[*].HGVS_P"   > dinaric.table
srun java -jar  /nbi/software/testing/snpEFF/4.3/x86_64/SnpSift.jar extractFields baltic.missense.vcf CHROM POS REF ALT AC AN  "ANN[*].HGVS_P"   > baltic.table
srun java -jar  /nbi/software/testing/snpEFF/4.3/x86_64/SnpSift.jar extractFields pannonean.missense.vcf CHROM POS REF ALT AC AN  "ANN[*].HGVS_P"   > pannonean.table
srun java -jar  /nbi/software/testing/snpEFF/4.3/x86_64/SnpSift.jar extractFields se_carpathian.missense.vcf CHROM POS REF ALT AC AN  "ANN[*].HGVS_P"   > se_carpathian.table
srun java -jar  /nbi/software/testing/snpEFF/4.3/x86_64/SnpSift.jar extractFields w_carpathian.missense.vcf CHROM POS REF ALT AC AN  "ANN[*].HGVS_P"   > w_carpathian.table

paste baltic.table w_carpathian.table se_carpathian.table dinaric.table pannonean.table > bal_wcarp_secarp_din_pan.table

awk '{print " "$1" "$2" "$3" "$4" "$5" "$6" "$12" "$13" "$19" "$20" "$26"  "$27" "$33" "$34" "$14}' bal_wcarp_secarp_din_pan.table > bal_wcarp_secarp_din_pan.format
