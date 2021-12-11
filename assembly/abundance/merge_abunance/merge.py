#usage:
#python merge.py input_filename outfile pathI

 
import sys,os
import pandas as pd
filename = sys.argv[1]
outfile = sys.argv[2]
pathI = sys.argv[3]
os.chdir(sys.argv[3])
units = pd.read_table(filename, dtype=str).set_index(["SRR"], drop=False)
input_list=[]
for i in units.itertuples():
	input_dir = "%s" % (i.SRR)
	input_list.append(input_dir)
#print input_list

counts = [pd.read_table(f, index_col=0, usecols=[0, 4], header=None, skiprows=0)
          for f in input_list]
#print counts
for t, (SRR) in zip(counts, units.index):
	t.columns = ["%s" % (SRR)]
matrix = pd.concat(counts, axis=1) 
matrix.index.name = "protein"

#print matrix.index
matrix.to_csv(outfile, sep="\t")
