Use a combination of bcftools and snpSIFT to collect data from the whole genome sequenced vcf into a format that can be passed onto FAAD.
 
source bcftools-1.2
  
## Extract the right scaffolds and contrasts and populations for your test
 
bcftools view   -S diploid_populations.txt  -r regions samples.ann.vcf.gz > regions_diploids.vcf
bcftools view   -S tetraploid_populations.txt -r regions samples.ann.vcf.gz > regions_tetraploids.vcf


----- 
 
source snpEFF-4.3
source jre-7.21

## Take the right fields from the annotated files
 
java -jar  SnpSift.jar extractFields regions_diploids.vcf CHROM POS REF ALT AC AN  "ANN[*].HGVS_P"   > diploid.table

java -jar  SnpSift.jar extractFields regions_tetraploids.vcf  CHROM POS REF ALT AC AN  "ANN[*].HGVS_P"  >  tetraploid.table

## Create a big table containing all the information

paste diploid.table tetraploid.table > diploid_tetraploid_contrast.table

## Manipulate this table to be readable by the python script

awk '{print " "$1" "$2" "$3" "$4" "$5" "$6" "$12" "$13" "$14}' diploid_tetraploid_contrast.table > diploid_tetraploid_contrast.txt
 
