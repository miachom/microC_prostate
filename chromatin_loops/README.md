This script detects chromatin loops and TAD borders using python package chromosight from the HiC matrices in cool format files. 
These commands will each generate 3 files in the results directory-

*loops.tsv (contains genomic coordinates, bin ids and correlation scores for the pattern identified

*loops.json (contains the windows of the same size as the kernel used around the patterns from pattern.txt)

*borders.json (contain images of the matrix regions around each detected border)

The command uses ```-m100000 -M500000000 -p0.35``` to detect all chromosome loops with sizes between 100kb and 50mb using 8 parallel threads.  
