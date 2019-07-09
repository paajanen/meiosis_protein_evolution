#!/usr/bin/env python
import re
import string
import sys
from itertools import combinations



inp=(sys.argv[1]) ## input file in FreqSum format

## This part reads the data into memory in a form that can be used easily for calculating Hamming distances.

with open ('{0}'.format(inp),'r') as  tets:
    lines_tet = tets.readlines()  
    gene_length=len(lines_tet)
    names = string.split(lines_tet[0])   ## takes the names of the samples
    haplotypes_d= [[] for _ in range(80)] ## list collecting the diploid haplotypes, we had 80 samples
    haplotypes_t = [[] for _ in range(40)] ## list collecting the tetraploid haplotypes, we had 40 samples
    for i in range(1,len(lines_tet)):
        fields=string.split(lines_tet[i])  
        for k in range(4,84): # the normal range of all diploids in our sampling
            haplotypes_d[k-4].append(fields[k])
        for j in range(84,124): ## the normal range of all tetraploids in our sampling
            haplotypes_t[j-84].append(fields[j])
    
## We calculate all the pairwise Hamming distances of the diploids, and take the maximum and call that a diameter.

    comb = list(combinations(range(0,80), 2)) ## enumerating the pairs of diploids
    hamming_distance=[[0] for _ in range(3160)] ## creating a list the length of the number of pairs
    for com in comb:
        ind=int(comb.index(com))
        for i in range(gene_length-1):
            if haplotypes_d[com[0]][i]=='0' and int(haplotypes_d[com[1]][i])==2 and haplotypes_d[com[0]][i]!='-1' and haplotypes_d[com[1]][i]!='-1':
                hamming_distance[ind][0]=hamming_distance[ind][0]+1
            if int(haplotypes_d[com[0]][i])==2 and int(haplotypes_d[com[1]][i])==0 and haplotypes_d[com[0]][i]!='-1' and haplotypes_d[com[1]][i]!='-1':
                hamming_distance[ind][0]=hamming_distance[ind][0]+1
    diameter=max(hamming_distance)
    pair=hamming_distance.index(diameter)
    local=min(hamming_distance)
    pair_local=hamming_distance.index(local)
    print max(hamming_distance), comb[pair], names[comb[pair][0]+4], names[comb[pair][1]+4], 'max diploids'
    #print min(hamming_distance), comb[pair_local], names[comb[pair_local][0]+4], names[comb[pair_local][1]+4],'min diploid'
    max_dip_1=names[comb[pair][0]+4]
    max_dip_2=names[comb[pair][1]+4]
    max_dip_1_index=comb[pair][0]+4
    max_dip_2_index=comb[pair][1]+4

## We calculate all the pairwise Hamming distances of the tetraploids, and take the maximum and call that a diameter.

    comb_t = list(combinations(range(0,40), 2))
    hamming_distance=[[0] for _ in range(780)]
    for com in comb_t:
        ind=int(comb_t.index(com))
        for i in range(gene_length-1):
            if haplotypes_t[com[0]][i]=='0' and int(haplotypes_t[com[1]][i])>2 and haplotypes_t[com[0]][i]!='-1' and haplotypes_t[com[1]][i]!='-1':
                hamming_distance[ind][0]=hamming_distance[ind][0]+1
            if int(haplotypes_t[com[0]][i])>2 and int(haplotypes_t[com[1]][i])==0 and haplotypes_t[com[0]][i]!='-1' and haplotypes_t[com[1]][i]!='-1':
                hamming_distance[ind][0]=hamming_distance[ind][0]+1
    diameter=max(hamming_distance)
    pair=hamming_distance.index(diameter)
    local=min(hamming_distance)
    pair_local=hamming_distance.index(local)
    print max(hamming_distance), comb_t[pair], names[comb_t[pair][0]+85], names[comb_t[pair][1]+85], 'max tetraploids'

    
## We calculate all the pairwise Hamming distances of diploids to tetraploids, and take the maximum and call that the maximum distance.
   
    comb_joint=[(d, t) for d in range(0,80) for t in range(0,39)]
    hamming_distance=[[0] for _ in range(3120)]
    for com in comb_joint:
        ind=int(comb_joint.index(com))
        for i in range(gene_length-1):
            if int(haplotypes_d[com[0]][i])==0 and int(haplotypes_t[com[1]][i])>2 and haplotypes_d[com[0]][i]!='-1' and haplotypes_t[com[1]][i]!='-1':
                hamming_distance[ind][0]=hamming_distance[ind][0]+1
            if int(haplotypes_d[com[0]][i])==2 and int(haplotypes_t[com[1]][i])==0 and haplotypes_d[com[0]][i]!='-1' and haplotypes_t[com[1]][i]!='-1':
                hamming_distance[ind][0]=hamming_distance[ind][0]+1
    
    diameter=max(hamming_distance)
    pair=hamming_distance.index(diameter)
    local=min(hamming_distance)
    pair_local=hamming_distance.index(local)
    print max(hamming_distance), comb_joint[pair], names[comb_joint[pair][0]+4], names[comb_joint[pair][1]+85], 'max mixed'
    #print min(hamming_distance), comb_joint[pair_local], names[comb_joint[pair_local][0]+4], names[comb_joint[pair_local][1]+85],'min mixed'
    max_mixed_d=names[comb_joint[pair][0]+4]
  
    
## If the maximum distance from diploid to tetraploid is the same as the one giving maximum diameter of the diploids, as sanity check we take the distance from the other part of the diploid pair to the tetraploid, to get a 'generic distance' in case the diploid was very diverged.    
    mixed_index=list(comb_joint).index((max_dip_2_index,comb_joint[pair][1]))
    if max_mixed_d==max_dip_1:
        mixed_index=list(comb_joint).index((max_dip_2_index,comb_joint[pair][1]))
        print hamming_distance[mixed_index], (max_dip_2_index-4,comb_joint[pair][1]), names[max_dip_2_index], names[comb_joint[pair][1]+85], 'general mixed'
    if max_mixed_d==max_dip_2:
        mixed_index=list(comb_joint).index((max_dip_1_index,comb_joint[pair][1]))
        print hamming_distance[mixed_index], (max_dip_1_index-4,comb_joint[pair][1]), names[max_dip_1_index], names[comb_joint[pair][1]+85], 'general mixed'
                
        
