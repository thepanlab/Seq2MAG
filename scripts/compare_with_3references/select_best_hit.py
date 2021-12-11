import os

dist = "specie883_Pasolli_et_al.Cell_0.8.dist"
with open(dist) as f1_obj:
    lines = f1_obj.readlines()
    
distance_dic = {}
hit_dic = {}

for line in lines:
    items = line.strip().split("\t")
    distance = items[2]
    query = items[1]
    reference = items[0]
    if distance_dic.get(query) == None:
        distance_dic[query] = distance
        hit_dic[query] = [reference, query]
    else:
        if distance < distance_dic[query]:
            distance_dic[query] = distance
            hit_dic[query] = [reference, query]
            
fout = open("species883_Pasolli_et_al_best_hit.ANI",'w')
for key, value in hit_dic.items():
    q = value[1]
    r = value[0]
    qname = os.path.basename(q)[:-3]
    rname = os.path.basename(r)[:-3]
    fout.write("ANIcalculator -genome1fna " + r + "  -genome2fna " + q + " -outfile " + qname+"_" + rname +".gANI" + " -outdir gANI_tmp"+"\n")

