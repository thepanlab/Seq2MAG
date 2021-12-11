# python script.py cluster_files output_file_name transpose_file_name
#import pickledb

#lengthDB = pickledb.load("protein_length.db",False)

#lengthDB.get(key)
import glob
import sys
import numpy

#kaver = 6724674.5

clusters_in = sys.argv[1]

output_file = sys.argv[2]
trans_file = sys.argv[3]

'''
clusters_in = "Cluster_113164_proteins"
output_file = "cluster113164_read_counts"
trans_file = "cluster113164_read_counts_transpose"
'''
def transpose(file_in,file_out):
    matrix = []
    with open(file_in) as ft_obj:
        for l in ft_obj.readlines():
            matrix.append(l.strip().split("\t"))
    ft_obj.close()
    trans = numpy.transpose(matrix)
    
    ft_out = open(file_out, 'w')
    for m in trans:
        ft_out.write("\t".join(map(str,m)) + "\n")
    
    ft_out.close()


def get_subjects(proteins_in_clus):
    protein_list = proteins_in_clus.split(",")
    subject_dic = {}
    for p in protein_list:
        infors = p.split("_")
        subject = infors[0]
        if subject in subject_dic.keys():
            subject_dic[subject].append(p)
        else:
            subject_dic[subject] = [p]
    return subject_dic

def length_dic(subjectID):
    protein_length = {}
    len_dir = "/work/TEDDY/cluster/scaffolds/737protein_length/"
    len_file = len_dir + subjectID + "_protein_length"
    with open(len_file) as f2_obj:
        for eachl in f2_obj.readlines():
            cols = eachl.strip().split("\t")
            protein_length[cols[0]] = cols[1]
    f2_obj.close()
    return protein_length

def get_sample_counts(sample_file_name,sid):
    counts_dic = {}
    with open(sample_file_name) as f4_obj:
        for reads_line in f4_obj.readlines():
            colls = reads_line.split("\t")
            scaffold_name = sid + "_" + colls[0]
            counts_dic[scaffold_name] = colls[2]
            #total = int(colls[1]) * 2 
    f4_obj.close()
    #counts_dic['total_counts'] = total  
    return counts_dic
        
## get subject list; total 737
subjects = []
with open("/work/TEDDY/737_subjects") as f3_obj:
    for x in f3_obj.readlines():
        subjects.append(x.strip())
f3_obj.close()

## get all clusters; key is the cluster id; value is the dictionary;
## the sub_dictionary is key with subject id, value with protein list in this subject
clusters = {}
with open(clusters_in) as f1_obj:
    clus_lines = f1_obj.readlines() # [1:3]
fout = open(output_file,'w')
fout.write("clusters\t") # because we will do transpose
for clus in clus_lines:
    # get cluster group ID and proteins in each cluster
    sp_clus = clus.strip().split("\t")
    clusID = sp_clus[0]
    proteins = sp_clus[1]
    proteins_each_subject = get_subjects(proteins)
    clusters[clusID] = proteins_each_subject
    #fout.write(clusID + "\t")
f1_obj.close() 
#fout.write("\n")
sort_keys = sorted(clusters.keys())
fout.write("\t".join(sort_keys) + "\n")

##  parse each subject
for a_subject in subjects:
    protein_len = length_dic(a_subject)
    counts_dir = "/work/TEDDY/abundance/finished/" + a_subject +"/*.rpkm"
    samples = glob.glob(counts_dir)
    #fout.write("clusterID" + "\t".join(samples)+ "\n")
    for s in samples:
        clus_counts = []
        sample_counts = get_sample_counts(s,a_subject)
        #ks = sample_counts['total_counts']
        #si = ks/kaver
            
        for cid in sort_keys:
            protein_dic = clusters[cid]
            #fout.write(cid + "\t")
            cid_sum = 0
            if a_subject in protein_dic.keys():
                proteins_used = protein_dic[a_subject]
                for apro in proteins_used:
                    scaff = "_".join(apro.split("_")[:-1])
                    scaff_len = float(apro.split("_")[4])
                    #protein_reads = int(sample_counts[scaff]) / scaff_len * int(protein_len[apro]) * int(ks) / float(kaver)
                    protein_reads = float(int(sample_counts[scaff]) / scaff_len * int(protein_len[apro]))

                    cid_sum = cid_sum + protein_reads
            else:
                cid_sum = 'NA'
            clus_counts.append(cid_sum)
        fout.write(s + "\t" + "\t".join(map(str,clus_counts)) + "\n")
fout.close()

transpose(output_file, trans_file)
    
