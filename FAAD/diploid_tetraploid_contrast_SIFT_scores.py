#!/usr/bin/env python
import re
import string
import sys
import gzip

Amino_acid_code={'Gly':'G','Ala':'A','Val':'V','Leu':'L','Ile':'I','Phe':'F','Trp':'W','Tyr':'Y','Asp':'D','Asn':'N','Glu':'E','Lys':'K','Gln':'Q','Met':'M','Ser':'S','Thr':'T','Cys':'C','Pro':'P','His':'H','Arg':'R'}

exponent=float((sys.argv[1]))
scaffold=(sys.argv[2]) ## and is the  infile to be used is from the prepared file from step 1, called diploid_tetraploid_contrast.txt 


with open ('/path_to_infile/'.format(scaffold),'r') as  Population_genomics:
    with open ('/path_to_SIFT_scores/{0}_SIFT_SCORES.txt'.format(scaffold),'r') as Annotated: 
        
## First we find a SIFT score for each locus.
        
        Loci=[]
        Annotations=[]
        for bline in Annotated:
            fields=string.split(string.strip(bline))
            Loci.append(fields[11])
            if fields[8]=='NONSYNONYMOUS':
                Annotations.append([fields[6],fields[8],fields[11],fields[12], fields[16]])
            if fields[9]=='NONSYNONYMOUS':
                Annotations.append([fields[6],fields[8],fields[12],fields[13], fields[17]])            
            if fields[10]=='NONSYNONYMOUS':
                Annotations.append([fields[6],fields[8],fields[13],fields[14], fields[18]])
        SIFT_dict=dict(zip(Loci,Annotations))

## Secondly we calculate the population genomics measures and combine it with the SIFT score.

        print "Chrom", "loci", "Amino Acid change", "Derived in diploid", "Derived in tetraploid", "DAP_DAF1", "DAP_DAF2", "SIFT score", "SIFT*DAP_DAF1", "SIFT*DAP_DAF2", "Absolute Difference", "Severity" 
        for aline in Pop1:
            fields=string.split(string.strip(aline))
            amino_acid=re.search(r'.{2}[\D]{3}(\d*)[\D]{3}',fields[8]) ## Recod that change in amino acids.
            if amino_acid:
                amino_acid_locus=amino_acid.group(1)
                derived=float(float(fields[4])+float(fields[6])) ## calculate how many derived alleles are there in total.
                if float(fields[5])>80 and float(fields[7])>80 and derived!=0:  ## Filter for half the depth, and for the fact that there are differences
                    dap=float((float(fields[4])/derived)**exponent+float(float(fields[6])/derived)**exponent) #Calculates the allele purity, here exponent present
                    daf_1=float(float(fields[4])/float(fields[5])) # derived allele frequency for diploids
                    daf_2=float(float(fields[6])/float(fields[7])) # deribed allele frequence for tetraploids
                    if amino_acid_locus in SIFT_dict.keys():  ## If this locus is found with a SIFT score, go and calculate.
                        if SIFT_dict[str(amino_acid_locus)][3][0]=='0' or SIFT_dict[str(amino_acid_locus)][3][0]=='1' :
                            SIFT=1-float(SIFT_dict[str(amino_acid_locus)][3])    ## Recalibrate the SIFT scale               
                            dap_daf1=float(dap*daf_1)
                            dap_daf2=float(dap*daf_2)
                            print fields[0],fields[1], fields[8], fields[4], fields[6], dap_daf1, dap_daf2, SIFT,float(dap_daf1*SIFT), float(dap_daf2*SIFT), abs(float(dap_daf1*SIFT)-float(dap_daf2*SIFT)), SIFT_dict[str(amino_acid_locus)][4]
                            
## Prints the scores.   
