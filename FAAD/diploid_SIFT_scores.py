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

        print "Chrom", "loci", "Amino Acid change","DAP_DAF1", "DAP_DAF2", "DAP_DAF3", "DAP_DAF4", "DAP_DAF5","SIFT score", "SIFT*DAP_DAF1", "SIFT*DAP_DAF2","SIFT*DAP_DAF3", "SIFT*DAP_DAF4","SIFT*DAP_DAF5", "Severity" 
        for aline in Pop1:
            fields=string.split(string.strip(aline))
            amino_acid=re.search(r'.{2}[\D]{3}(\d*)[\D]{3}',fields[14])
            if amino_acid:
                amino_acid_locus=amino_acid.group(1)
                derived=float(float(fields[4])+float(fields[6])+float(fields[8])+float(fields[10])+float(fields[12])) # check that the site is non-redundant
                if float(fields[5])>16 and float(fields[7])>16 and float(fields[9])>16 and float(fields[11])>16 and float(fields[13])>16 and derived!=0: # filter sites that are bad quality
                    dap=float((float(fields[4])/derived)**exponent+float(float(fields[6])/derived)**exponent+float(float(fields[8])/derived)**exponent+float(float(fields[10])/derived)**exponent+float(float(fields[12])/derived)**exponent) #Calculates the allele purity, here exponent present
                    daf_1=float(float(fields[4])/float(fields[5])) # derived allele frequency for diploid population 1
                    daf_2=float(float(fields[6])/float(fields[7])) # derived allele frequency for diploid population 2
                    daf_3=float(float(fields[8])/float(fields[9])) # derived allele frequency for diploid population 3
                    daf_4=float(float(fields[10])/float(fields[11])) # derived allele frequency for diploid population 4
                    daf_5=float(float(fields[12])/float(fields[13])) # derived allele frequency for diploid populatoin 5
                    if amino_acid_locus in SIFT_dict.keys():
                        if SIFT_dict[amino_acid_locus][3][0]=='0' or SIFT_dict[amino_acid_locus][3][0]=='1' :
                            SIFT=1-float(SIFT_dict[amino_acid_locus][3])                   
                            dap_daf1=float(dap*daf_1)
                            dap_daf2=float(dap*daf_2)
                            dap_daf3=float(dap*daf_3)
                            dap_daf4=float(dap*daf_4)
                            dap_daf5=float(dap*daf_5)
                            print fields[0],fields[1], fields[14], dap_daf1, dap_daf2, dap_daf3, dap_daf4, dap_daf5, SIFT ,float(dap_daf1*SIFT), float(dap_daf2*SIFT), float(dap_daf3*SIFT),float(dap_daf4*SIFT),float(dap_daf5*SIFT), SIFT_dict[amino_acid_locus][4]
