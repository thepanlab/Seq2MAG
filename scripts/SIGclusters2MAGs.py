#python
import pickledb
import json
import math

# step 1 creat protein_class pickleDB
cluster_file = "protein_family_filename"
db_protein_class = "protein_class_core_peri.db"
db_protein_cluster = "protein_clusterID.db"

protein_class = {}
protein_cluster = {}
with open(cluster_file) as cluster_obj:
    for line in cluster_obj.readlines()[1:]:
        groups = line.strip().split("\t")
        subject_number = int(groups[1])
        clusterID = "_".join(groups[0].split())[1:]
        proteins = groups[3].split(",")
        for protein in proteins:
            protein_cluster[protein] = clusterID
            if subject_number > 737 * 0.5:
                protein_class[protein] = "core"
            else:
                protein_class[protein] = "peri"


def creat_pickleDB (dic,file_out):
    n = len(dic)
    x =0 
    with open(file_out,'w') as out_obj:
        out_obj.write('{')
        for key, value in dic.items():
            x += 1
            if x != n:
                out_obj.write('"' + key + '": "' + value + '", ' )
            elif x == n:
                out_obj.write('"' + key + '": "' + value + '"')
        out_obj.write('}')
        
        
creat_pickleDB(protein_class, db_protein_class)
creat_pickleDB(protein_cluster, db_protein_cluster)       

## step2 creat scaffold --> proteins json database

#subject_proteins = open('/work/TEDDY/cluster/scaffolds/737protein_length/998650_protein_length')
subject_proteins = open("737_subject_proteins")
#subject = "998650"
proteins_in_scaffold = {}
for each_line in subject_proteins.readlines():
    items = each_line.strip().split("\t")
    protein_tmp = items[0]
    scaffold = '_'.join(protein_tmp.split("_")[:-1])
    #print (scaffold)
    length = int(items[1]) # here is the cds length, so the length > 50 is not right.
    # the protein length is >50 aa
    if length > 149:
#        print (proteins_in_scaffold.get(scaffold))
        if proteins_in_scaffold.get(scaffold) == None:
            proteins_in_scaffold[scaffold] = [protein_tmp]
            #print(proteins_in_scaffold[scaffold])
        elif proteins_in_scaffold.get(scaffold) != 'None':
            proteins_in_scaffold[scaffold].append(protein_tmp)
import json
with open("all_scaffold2proteins.json",'w') as fp:
    json.dump(proteins_in_scaffold,fp,indent=2)
        

# step3 creat MAGs--> scaffold json database

all_MAGs = open("MAGs2scaffold")
scaffold_in_MAGs = {}
for each in all_MAGs.readlines():
    cols = each.strip().split("\t")
    MAG = ".".join(cols[0].split(".")[:-1])
    subject = MAG.split("_")[0]
    scaffold_tmp = subject + "_" + cols[1]
    if scaffold_in_MAGs.get(MAG) == None:
        scaffold_in_MAGs[MAG] = [scaffold_tmp]
    elif scaffold_in_MAGs.get(MAG) != None:
        scaffold_in_MAGs[MAG].append(scaffold_tmp)
import json
with open("all_MAGs2scaffold.json",'w') as fp2:
    json.dump(scaffold_in_MAGs,fp2,indent=2)

# step4 creat cluster--> q-value/coeff/p-value json database

glmm = open("/work/TEDDY/MAG_core_protein/glmm_final/pre_after_IA_result.tsv")
glmm_result = {}
for each_row in glmm.readlines():
    each_cols = each_row.strip().split("\t")
    clusID = each_cols[0]
    glmm_result[clusID] = [each_cols[1],each_cols[2],each_cols[3]]
with open("glmm_result.json",'w') as js:
    json.dump(glmm_result,js,indent=2)
    
    
#load databases

db_class = pickledb.load(db_protein_class,False)
db_clusterID = pickledb.load(db_protein_cluster,False)

with open('all_MAGs2scaffold.json') as fd1:
    db_MAGs = json.load(fd1)

with open('all_scaffold2proteins.json') as fd2:
    db_scaffold = json.load(fd2)

with open('glmm_result.json') as fd3:
    db_glmm = json.load(fd3)
fd1.close()
fd2.close()
fd3.close()


MAG_sum_out = open("MAG_sum_out_E-6_coeff_Ncluster.tsv",'w')
MAG_sum_out.write("MAG_name\tall_numbers\tmissing\tcore_Pnumbers\tPup_core\tPdown_core\tPno_sig_core\tcore_Cnumbers\tCup_core\tCdown_core\tCno_sig_core\n")

def get_MAG_info (mag_name):
    file_out = open("/work/TEDDY/MAG_core_protein/high_mags_tmp/" + mag_name + ".proteins.tsv",'w')
    scaffolds_in_this_MAG = db_MAGs[mag_name]
    core_numbers = 0
    peri_numbers = 0
    all_numbers = 0
    up = 0
    down = 0
    not_sig = 0
    core_up_cluster =[]
    core_down_cluster = []
    core_not_sig_cluster = []
    core_cluster = []
    for s in scaffolds_in_this_MAG:
        if db_scaffold.get(s) != None:
            proteins_in_this_scaffold = db_scaffold[s]
            
            for p in proteins_in_this_scaffold:
                all_numbers += 1
                core_peri_info = db_class.get(p)
                p_clusterID = db_clusterID.get(p)
                #print(mag_name,p,core_peri_info,p_clusterID)
                #file_out.write(mag_name + "\t" + p + "\t" + str(core_peri_info) + "\t" + str(p_clusterID) + "\n")
                if core_peri_info == 'core':
                    core_numbers += 1
                    core_cluster.append(p_clusterID)
                    pvalue = float(db_glmm[p_clusterID][0])
                    qvalue = float(db_glmm[p_clusterID][1])
                    coeff = float(db_glmm[p_clusterID][2])
                    # parse up and down clusters
                    
                    # condition q < E-6
                    if qvalue < 0.000001 and coeff > 0:
                        up += 1
                        core_up_cluster.append(p_clusterID)
                    elif qvalue < 0.000001 and coeff <0:
                        down += 1
                        core_down_cluster.append(p_clusterID)
                    elif qvalue >= 0.000001:
                        not_sig += 1
                        core_not_sig_cluster.append(p_clusterID)
                    '''
                    
                    # condition q < E-4 and |coeff| > log2(1.5)
                    if qvalue < 0.0001 and coeff > math.log2(1.5):
                        up += 1
                    elif qvalue < 0.0001 and coeff < -(math.log2(1.5)):
                        down += 1
                    else:
                        not_sig += 1
                    '''
                    file_out.write(mag_name + "\t" + p + "\t" + str(core_peri_info) + "\t" + str(p_clusterID) + "\t" + str(pvalue) + "\t" + str(coeff) + "\n")
                elif core_peri_info == 'peri':
                    peri_numbers += 1
    MAG_sum_out.write(mag_name + "\t" + str(all_numbers) + "\t" + str(all_numbers-core_numbers-peri_numbers) + "\t" + str(core_numbers) + "\t")
    MAG_sum_out.write(str(up) + "\t" + str(down) + "\t" + str(not_sig) + "\t")
    MAG_sum_out.write(str(len(set(core_cluster)))+ "\t" + str(len(set(core_up_cluster))) + "\t" +  str(len(set(core_down_cluster))) + "\t" + str(len(set(core_not_sig_cluster))) + "\n")
    file_out.close()
#get_MAG_info('200318_bin.103')
#c = 0
with open('high_quality_final') as mag_file:
    for aline in mag_file.readlines():
        #c += 1
        each_mag_name = aline.split()[0]
        #if isinstance(c/100,int): 
        #    print(c)
        get_MAG_info(each_mag_name)
    
mag_file.close()
MAG_sum_out.close()
