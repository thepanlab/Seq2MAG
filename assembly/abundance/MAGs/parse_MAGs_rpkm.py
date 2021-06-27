#usage 'python parse_MAGs_rpkm.py count_result total_reads_file srrID subject_mags_tmp'

import re
import sys

count_result = sys.argv[1]
total_reads_file = sys.argv[2]
srr = sys.argv[3]
mags = sys.argv[4]
outfile = mags + "_" + srr + "_rpkm.tsv"

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


scaffold_length = {}
scaffold_counts = {}

for row in rows:
	groups = row.strip().split()
	node = groups[0].split('_')
	counts = int(groups[2])
	length = int(node[3])
	scaffold_length[groups[0]] = length
	scaffold_counts[groups[0]] = counts

file_object.close()

with open(mags) as mags_object:
    bins = mags_object.readlines()

f = open (outfile,'w')

totals = 2*int(total(srr))
bin_length = {}
bin_counts = {}

for bin in bins:
	cols = bin.strip().split()
	scaffold = cols[1]
	bin_name = cols[0]
	if bin_name not in bin_length.keys():
		bin_length[bin_name] = scaffold_length[scaffold]
		bin_counts[bin_name] = scaffold_counts[scaffold]
	else:
		bin_length[bin_name] = bin_length[bin_name] + scaffold_length[scaffold]
		bin_counts[bin_name] = bin_counts[bin_name] + scaffold_counts[scaffold]

		
for key, value in bin_length.items():
	rpkm = bin_counts[key] * 1e9/(totals*value)
	f.write(key + "\t" + str(bin_counts[key]) + "\t" + str(value) + "\t" + str(totals) + "\t" + str(rpkm) + "\n")	

"""
with open(outfile,'w') as out_object:
    for row in rows:
        groups = row.strip().split()
        node = groups[0].split('_')
        counts = int(groups[2])
        length = int(node[3])
        totals = int(total(srr))
        rpkm = counts/(length/1000*totals/1e6)
        
        out_object.write(groups[0] + "\t" + str(totals) + "\t" + str(counts) + "\t" + str(length) + "\t"+ str(rpkm) + "\n")
        
"""
