import sys,os
import re	
import glob

#usage 'python parse_abundance_rpkm.py cluster_input rpkm_output'

F_in = sys.argv[1]
f4 = open(sys.argv[2],'w')


with open(F_in) as file_object:
	lines = file_object.readlines()

def my_grep(Aprotein):
	id = Aprotein[0:7]
	file_list = glob.glob('/work/TEDDY/abundance/rpkm/'+str(id)+'*pkm.tsv')
	if len(file_list) == 1:
		file = file_list[0]
		#print file
		with open(file) as f2:
			rpkm_line = [next(f2)]
			for line in f2:
				contents = line.strip().split("\t")
				if contents[0] in Aprotein:
					rpkm_line.append(line)
	else:
		print "no " + id
		rpkm_line=[]
	
	return rpkm_line	
def sum_each(subject,F,subject_list):
	count = {}
	n = 0
	final_subject = re.sub("_","",subject)
	titles = re.sub("\n","",subject_list[0][0])
	ngroups = titles.strip().split("\t")
	for i in range(1,len(ngroups)):
		count[str(i)] = [ngroups[i]]
	for each_protein in subject_list:
		rm = re.sub("\n","",each_protein[1])
		ortholog = rm.strip().split("\t")
		for m in range(1,len(ortholog)):
			count[str(m)].append(ortholog[m])
	f4.write(F + "\t" + final_subject + "\t")
	for val in count.values():
		month = re.sub("_","",val[0][0:2])
		srr = re.sub("_","",val[0][-16:-4])#new
		srr = re.sub("\.","",srr)#new
		rpkm = map(float, val[1:])
		#print rpkm
		#f4.write("(" + str(month) + ","+ str(sum(rpkm)) + ");")
		f4.write("(" + str(month) + "," + str(srr) + "," + str(sum(rpkm)) + ");")
		rpkm = ""

	f4.write("\n")
			

for line in lines:
	each_subject = {}
	groups = line.strip().split("\t")
	F = groups[0]
	proteins = groups[1].split(",")
	for protein in proteins:
		key = protein[0:7]
		new_list = my_grep(protein)
		if key in each_subject.keys() and len(new_list) >= 1:
			each_subject[key].append(new_list)
		elif key not in each_subject.keys() and len(new_list) >= 1:
			each_subject[key] = [new_list]
	for i in each_subject.keys():
		#print i, F, each_subject[i]
		sum_each(i,F,each_subject[i])

