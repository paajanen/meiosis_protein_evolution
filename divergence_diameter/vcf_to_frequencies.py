#!/usr/bin/env python
import re
import string
import sys
import gzip

for aline in sys.stdin:
    if aline[0:2]=='##': #Read the header of the vcf, but don't print it
        continue
    if aline[0:2]=='#C': #Take the fields line, modify and print it in the required format
        ref_fields=string.split(string.strip(aline))
        no_of_samples=len(ref_fields)-9
        header=[ref_fields[0],ref_fields[1],ref_fields[3],ref_fields[4]]
        for i in range(9,9+no_of_samples):
            header.append(ref_fields[i]+"("+str(4)+")")
        print "\t".join(map(str, header))
    else: # Parse the vcf and print it in required format based on genotypes
        vcf_fields=string.split(string.strip(aline))
        alt_counts=[vcf_fields[0][9],vcf_fields[1],vcf_fields[3],vcf_fields[4]]
        for i in range(9,9+no_of_samples):
            genotype=re.search(r'^([01/.]+):',vcf_fields[i])
            if genotype.group(1)=='0/0':
                alt_count=0
            if genotype.group(1)=='0/1':
                alt_count=1
            if genotype.group(1)=='1/1':
                alt_count=2
            if genotype.group(1)=='./.':
                alt_count=-1            
            if genotype.group(1)=='0/0/0/0':
                alt_count=0
            if genotype.group(1)=='0/0/0/1':
                alt_count=1
            if genotype.group(1)=='0/0/1/1':
                alt_count=2
            if genotype.group(1)=='0/1/1/1':
                alt_count=3
            if genotype.group(1)=='1/1/1/1':  
                alt_count=4          
            if genotype.group(1)=='./././.':
                alt_count=-1            
            alt_counts.append(alt_count)
        print "\t".join(map(str, alt_counts))
