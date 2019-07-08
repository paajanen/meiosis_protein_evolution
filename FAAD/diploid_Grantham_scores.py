## uses python-2.7.10 and biopython-1.69

#!/usr/bin/env python
import re
import string
import sys
import gzip
import Bio

from Bio.SubsMat import MatrixInfo    
blosum = MatrixInfo.blosum100
grant = MatrixInfo.grant


Amino_acid_code={'Gly':'G','Ala':'A','Val':'V','Leu':'L','Ile':'I','Phe':'F','Trp':'W','Tyr':'Y','Asp':'D','Asn':'N','Glu':'E','Lys':'K','Gln':'Q','Met':'M','Ser':'S','Thr':'T','Cys':'C','Pro':'P','His':'H','Arg':'R'}
Amino_acid_code_grant={'Gly':'Y','Ala':'L','Val':'V','Leu':'U','Ile':'W','Phe':'F','Trp':'T','Tyr':'O','Asp':'A','Asn':'N','Glu':'G','Lys':'I','Gln':'Q','Met':'M','Ser':'S','Thr':'E','Cys':'C','Pro':'P','His':'H','Arg':'R'}


exponent=float((sys.argv[1]))
scaffold=(sys.argv[2])


with open ('/path_to_infile/{0}'.format(scaffold),'r') as  Pop1:
    
    print "Chrom", "loci", "Amino Acid change", "DAP_DAF1", "DAP_DAF2", "DAP_DAF3","DAP_DAF4","DAP_DAF5",  "Grantham matrix score", "Grantham*DAP_DAF1", "Grantham*DAP_DAF2", "Grantham*DAP_DAF3","Grantham*DAP_DAF4","Grantham*DAP_DAF5"
    for aline in Pop1:
        fields=string.split(string.strip(aline))
        if len(fields)<15:
            continue
        if fields[14][0]!="p":
            continue
        if fields[14][0:5]=="p.Ter":
            continue
        else:
            derived=float(float(fields[4])+float(fields[6])+float(fields[8])+float(fields[10])+float(fields[12])) # check that the site is non-redundant
            if float(fields[5])>16 and float(fields[7])>16 and float(fields[9])>16 and float(fields[11])>16 and float(fields[13])>16 and derived!=0: # filter sites that are bad quality
                dap=float((float(fields[4])/derived)**exponent+float(float(fields[6])/derived)**exponent+float(float(fields[8])/derived)**exponent+float(float(fields[10])/derived)**exponent+float(float(fields[12])/derived)**exponent)
                daf_1=float(float(fields[4])/float(fields[5]))   # derived allele frequency for diploid population 1  
                daf_2=float(float(fields[6])/float(fields[7]))   # derived allele frequency for diploid population 2
                daf_3=float(float(fields[8])/float(fields[9]))   # derived allele frequency for diploid population 3
                daf_4=float(float(fields[10])/float(fields[11])) # derived allele frequency for diploid population 4
                daf_5=float(float(fields[12])/float(fields[13])) # derived allele frequency for diploid populatoin 5
                amino_acid_change=re.search(r'.{2}([\D]{3})\d*([\D]{3})',fields[14])
                if amino_acid_change:
                    if amino_acid_change.group(1)!=amino_acid_change.group(2):
                        pair_grant = (Amino_acid_code_grant[str(amino_acid_change.group(1))],Amino_acid_code_grant[str(amino_acid_change.group(2))])
                        if pair_grant not in grant:
                            grant_matrix_score= grant[(tuple(reversed(pair_grant)))]
                        else:
                            grant_matrix_score = grant[pair_grant] 
                        dap_daf1=float(dap*daf_1)
                        dap_daf2=float(dap*daf_2)
                        dap_daf3=float(dap*daf_3)
                        dap_daf4=float(dap*daf_4)
                        dap_daf5=float(dap*daf_5)
                        print fields[0],fields[1], fields[14], dap_daf1, dap_daf2, dap_daf3, dap_daf4, dap_daf5, grant_matrix_score,float(dap_daf1*grant_matrix_score), float(dap_daf2*grant_matrix_score), float(dap_daf3*grant_matrix_score),float(dap_daf4*grant_matrix_score),float(dap_daf5*grant_matrix_score)
                else: continue 
