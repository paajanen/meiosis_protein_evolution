#!/bin/bash -e
source python-2.7.10
source biopython-1.69
## Run the Python script with various exponent
 
diploid_tetraploid_contrast_Grantham_scores.py 2 diploid_tetraploid_contrast.txt  > dip_tet_2.scores
diploid_tetraploid_contrast_Grantham_scores.py 2.5 diploid_tetraploid_contrast.txt  > dip_tet_2.5.scores
diploid_tetraploid_contrast_Grantham_scores.py 3 diploid_tetraploid_contrast.txt  > dip_tet_3.scores
diploid_tetraploid_contrast_Grantham_scores.py 3.5 diploid_tetraploid_contrast.txt  > dip_tet_3.5.scores
