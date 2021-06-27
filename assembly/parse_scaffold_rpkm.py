#usage 'python parse_rpkm.py count_result total_reads_file srrID'

import re
import sys

count_result = sys.argv[1]
total_reads_file = sys.argv[2]
srr = sys.argv[3]
outfile = srr + "rpkm.txt"

def total(srr):
    with open(total_reads_file) as f:
        lines = f.readlines()

    for line in lines:
        if line.find(srr) != -1:
            reads = line.strip().split("\t")
            return reads[0]
    f.close()

with open(count_result) as file_object:
    rows = file_object.readlines()

with open(outfile,'w') as out_object:
    for row in rows:
        groups = row.strip().split()
        node = groups[0].split('_')
        counts = int(groups[2])
        length = int(node[3])
        totals = int(total(srr))
        rpkm = counts/(length/1000*totals/1e6)
        
        out_object.write(groups[0] + "\t" + str(totals) + "\t" + str(counts) + "\t" + str(length) + "\t"+ str(rpkm) + "\n")
        

