From the dip_tet_2.scores result sort the column number 11, assuming that you have a contrast of two things, tail refers to 1% outliers (check by wc -l the number of snps you get, and take 1%)
 
sort -k11 -n dip_tet_2.scores | tail -21085 | sort -k 1,1 -k 2,2n >  dip_tet_outliers_2_0.01.sorted
