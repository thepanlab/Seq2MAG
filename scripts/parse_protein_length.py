# python parse_protein_length.py gff_file output

import sys
gff = sys.argv[1]
#out = sys.argv[2]
#gff ="300220_scaffolds_min1000_predicted.gff"
subject = gff.split("_")[0]
out = subject + "_protein_length"
fout = open(out,'w')

with open(gff) as f1_obj:
	lines = f1_obj.readlines()

for line in lines:
	if 'NODE_' in line and 'CDS' in line:
		infors = line.strip().split("\t")
		scaffold = infors[0]
		id_all = infors[8].split(";")[0]
		id_number = id_all.split("_")[1]
		protein = subject + "_" + scaffold + "_" + id_number
		start = int(infors[3])
		end = int(infors[4])
		length = abs(start-end)
		fout.write(protein + "\t" + str(length) + "\n")