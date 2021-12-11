import sys

#usage 'python creat_RXN_dictionary.py input.blst out_reaction_dictionary.txt'

blst = sys.argv[1]
#blst = "besthit_scaffolds_min1000_protein.blst"
outfile = sys.argv[2]
#outfile = "test_reaction_dic.txt"
fasta = "/work/omicsbio/lizhang12/tools/metacyc/database/uniprot-seq-ids.fasta"

#remove duplications in a list
def remove_dup(list_x):
    return list(dict.fromkeys(list_x))

# reaction and EC dictionary; the key is protein name and the vaule is the reactions and ECs
reactions = {}
EC ={}

with open(fasta) as fasta_object:
    for row in fasta_object:
        if '>' in row:
            titles = row.strip().split(" ")
            protein = titles[0][1:]
            rxn = titles[1] # separate with ','
            reactions[protein] = rxn
            if len(titles) == 3:
                ec_number = titles[2]
                EC[protein] = ec_number
            elif len(titles) == 2 and 'RXN' in row:
                EC[protein] = 'NA'
fasta_object.close()

#print (reactions["tr|Q9WYS2|Oxidoreductase,_short_chain_dehydrogenase/reductase_family"])
#print (EC["tr|Q9WYS2|Oxidoreductase,_short_chain_dehydrogenase/reductase_family"])

select_query = {}
with open(blst) as blast_object:
    for line in blast_object:
        cols = line.strip().split("\t")
        query = cols[0]
        subject = cols[1]
        identity = float(cols[2])
        evalue = float(cols[10])
        
        #qcovs = cols[12]
        if identity >= 20 and evalue <= 0.0001:
            if query not in select_query.keys():
                #print (subject)
                select_query[query] = [reactions[subject]]
            #else:  #commit here if just select the best one
            #    select_query[query].append(reactions[subject])

fout = open(outfile,'w')

for p in select_query.keys():
    reactions_each_query = []
    for r_items in select_query[p]:
        all_rxn_in_value = r_items.split(",")
        reactions_each_query.extend(all_rxn_in_value)
    select_query[p] = remove_dup(reactions_each_query)
    fout.write(p + "\t" + ",".join(select_query[p]) + "\n")
