#
import os
import glob
path = '/work/TEDDY/binning/drep/compare/representive'
HGM_hit = "HGM_best_hit.gANI"
IGG_hit = "IGG_best_hit.gANI"
CELL_hit = "CELL_best_hit.gANI"

# key is representive name; value is a list with 2 items (1,HGM hit; 2, IGG hit)
dic_hit = {}
dic_hit_genome = {}

# get reprenstive names

for filename_tmp in glob.glob(os.path.join(path,'*.fa')):
    filename = os.path.basename(filename_tmp)
    dic_hit[filename] = ["NO","NO","NO"]
    dic_hit_genome[filename] = ["NO","NO","NO"]

# read HGM hit

def select_hit(row):
    items = row.strip().split("\t")
    genome1 = items[0]
    genome2 = items[1]
    ani1 = float(items[2])
    ani2 = float(items[3])
    af1 = float(items[4])
    af2 = float(items[5])
    if (ani1 > 95 or ani2 >95) and (af1 > 0.3 or af2 > 0.3):
        return genome1,genome2,"YES"
    else:
        return genome1,genome2,"NO"

count_HGM= 0
with open(HGM_hit) as f1_obj:
    for line in f1_obj.readlines()[1:]:
        referGenome, genomeID,hit_result = select_hit(line)
        if hit_result == "YES":
            count_HGM += 1
            dic_hit_genome[genomeID][0] = referGenome
        dic_hit[genomeID][0] = hit_result
        
    
count_IGG= 0
with open(IGG_hit) as f2_obj:
    for iline in f2_obj.readlines()[1:]:
        referGenome2,genomeID2,hit_result2 = select_hit(iline)
        if hit_result2 == "YES":
            count_IGG += 1
            dic_hit_genome[genomeID2][1] = referGenome2
        dic_hit[genomeID2][1] = hit_result2
        
count_CELL= 0
with open(CELL_hit) as f3_obj:
    for mline in f3_obj.readlines()[1:]:
        referGenome3,genomeID3,hit_result3 = select_hit(mline)
        if hit_result3 == "YES":
            count_CELL += 1
            dic_hit_genome[genomeID3][2] = referGenome3
        dic_hit[genomeID3][2] = hit_result3
print("HGM ",count_HGM)
print("IGG ", count_IGG)
print("CELL",count_CELL)

fout =open("883species_hit_HGM_IGG_CELL.0809.tsv",'w')
fout.write("species\tnewSpecies\tHGM\tIGG\tCELL\tHGM_hit\tIGG_hit\tCELL_hit\n")

n_all = 0
for g, hit_values in dic_hit.items():
    fout.write(g + "\t")
    if hit_values[0] == 'YES' or hit_values[1] == 'YES'or hit_values[2] == 'YES':
        n_all = n_all + 1
        fout.write("NO" + "\t")
    else:
        fout.write("YES" + "\t")6
    fout.write("\t".join(hit_values) + "\t" + "\t".join(dic_hit_genome[g]) + "\n")
    

print("all ", n_all)
        
    
    
    
    
    
    
    
    