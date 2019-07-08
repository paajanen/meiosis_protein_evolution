#### requires python-2.7.10 and biopython-1.69
 
#!/usr/bin/env python
import re
import string
import sys
import gzip
import Bio
 
from Bio.SubsMat import MatrixInfo   
blosum = MatrixInfo.blosum100
grant = MatrixInfo.grant
#blosum['D','N']
 
Amino_acid_code={'Gly':'G','Ala':'A','Val':'V','Leu':'L','Ile':'I','Phe':'F','Trp':'W','Tyr':'Y','Asp':'D','Asn':'N','Glu':'E','Lys':'K','Gln':'Q','Met':'M','Ser':'S','Thr':'T','Cys':'C','Pro':'P','His':'H','Arg':'R'}
Amino_acid_code_grant={'Gly':'Y','Ala':'L','Val':'V','Leu':'U','Ile':'W','Phe':'F','Trp':'T','Tyr':'O','Asp':'A','Asn':'N','Glu':'G','Lys':'I','Gln':'Q','Met':'M','Ser':'S','Thr':'E','Cys':'C','Pro':'P','His':'H','Arg':'R'}
 
 
exponent=float((sys.argv[1]))
scaffold=(sys.argv[2]) ## and is the  infile to be used is from the prepared file from step 1, called diploid_tetraploid_contrast.txt 
 
with open ('path_to_infile/'.format(scaffold),'r') as  Pop1:  

    print "Chrom", "loci", "AA", "AF/2x", "AF/4X",  "DAP_DAF1", "DAP_DAF2", "grant_matrix_score", "grant1", "grant2", "difference in scores" # Print the header line
    for aline in Pop1:
        fields=string.split(string.strip(aline))
        if fields[8][0]!="p":   ##check that you are having real amino acid changes
            continue
        else:
            derived=float(float(fields[4])+float(fields[6]))
            if float(fields[5])>80 and float(fields[7])>80 and derived!=0:  ## Filter for half the depth, and for the fact that there are differences
                dap=float((float(fields[4])/derived)**exponent+float(float(fields[6])/derived)**exponent) #Calculates the allele purity, here exponent present
                daf_1=float(float(fields[4])/float(fields[5])) # derived allele frequency for diploids
                daf_2=float(float(fields[6])/float(fields[7])) # deribed allele frequence for tetraploids
                amino_acid_change=re.search(r'.{2}([\D]{3})\d*([\D]{3})',fields[8]) ## Fine the change in amino acids
                if amino_acid_change:
                    pair_grant = (Amino_acid_code_grant[str(amino_acid_change.group(1))],Amino_acid_code_grant[str(amino_acid_change.group(2))]) ##turn these to one letter coder
                    if pair_grant not in grant:
                        grant_matrix_score= grant[(tuple(reversed(pair_grant)))] ## read the code from the grantham matrix info
                    else:
                        grant_matrix_score = grant[pair_grant]   ## and it depends on the order, as it upper triangular matrix
                else: continue                 
                dap_daf1=float(dap*daf_1)
                dap_daf2=float(dap*daf_2)
                print fields[0],"\t",fields[1],"\t", fields[8],"\t", fields[4],"\t", fields[6],"\t", dap_daf1,"\t", dap_daf2, "\t",grant_matrix_score,"\t",float(dap_daf1*grant_matrix_score),"\t", float(dap_daf2*grant_matrix_score), "\t",abs(float(dap_daf1*grant_matrix_score)-float(dap_daf2*grant_matrix_score))
 
##Calculates scores and prints the results line by line
 
