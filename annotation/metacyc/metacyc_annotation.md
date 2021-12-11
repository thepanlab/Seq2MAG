# MetaCyc Reaction annotation

```
A common early step in performing pathway analysis of genomes and metagenomes is to associate protein sequences
to MetaCyc reactions. The Pathway Tools software[1] infers such associations by using EC numbers, enzyme names, 
and Gene Ontology terms within protein annotations. Such annotations might be inferred using a variety of
sequence-analysis methods.

To aid researchers in associating sequences to MetaCyc reactions, MetaCyc enzymes that have a link to 
UniProt contain protein sequence information. It is possible to perform BLAST searches against MetaCyc proteins
with sequence information using the BLAST Search. 
In addition, each release of MetaCyc includes **a file that associates MetaCyc reaction IDs with the UniProt identifiers 
of enzymes known to catalyze those reactions. Note that not all MetaCyc reactions have EC numbers (because not all enzymes
have yet been assigned EC numbers), therefore EC numbers are not a comprehensive mechanism for associating sequences to 
reactions.** The file is called **uniprot‑seq‑ids.dat** and included in the MetaCyc data file distribution.[2,3,4]

references:
1.Karp P.D., et al. Pathway Tools version 19.0: Integrated Software for Pathway/Genome Informatics and Systems Biology
Briefings in Bioinformatics doi: 10.1093/bib/bbv079. (2015)
2.Caspi R., et al., The MetaCyc database of metabolic pathways and enzymes Nucleic Acids Research 46(D1):D633-9.(2018)
3.https://metacyc.org/
4.https://metacyc.org/MetaCycUserGuide.shtml; section "7  BLASTing Against MetaCyc"
```
## prepare database
We will use three database files in this github **(created May, 2020)**:

1.`uniprot-seq-ids.fasta`, this fasta file includes protein sequences with uniprot ID.

```
>sp|A5D7R8|Thymidine_kinase,_cytosolic THYKI-RXN,DURIDKI-RXN 2.7.1.21
MSCINLPNVLPGSPSKTRGQIQVILGPMFSGKSTELMRRVRRFQVAQYKCLVIKYAKDTRYSSLFSTHDRNTMEALPACLLRDVIQDAQRVAVIGIDEGQFFPDIVEFCENMANSGK
TVIVAALDGTFQRKAFGTILNLVPLAESVVKLTAVCMECFREAAYTKRLGVEKEVEVIGGADKYHSVCRLCYFKKASGQPAVLDSEENKENCPMTLGKPAEAPGVRKLFATHQIWQC
SQAN
```
2.`uniprot-seq-ids.dat`, this table link uniprot ID with Metacyc Reaction IDs.
```
 example:
 (GLYCEROL-DEHYDROGENASE-ACCEPTOR-RXN 1.1.99.22 "Q70JN9" "Q8L1D5")

The above 2 files are copied from Pathway tools. if you want to know how to get the uniprot-seq-ids.dat and 
uniprot-seq-ids.fasta, please see details in the following links:
http://bioinformatics.ai.sri.com/ptools/flatfile-format.html

```
3.`All_pathways_of_MetaCyc.txt`, this table link Reaction IDs wiht pathway IDs.

```
login to Metacyc database, and then creat your own specical smartable of 'All pathways of MetaCyc' (you can add columns 
in this samrtable to get more information of each pathway and download this smartable). 
https://metacyc.org/;
https://metacyc.org/PToolsWebsiteHowto.shtml#node_sec_6
```
### step1 Using DIAMOND to perform blast searches:
using 'uniprot-seq-ids.fasta' as blast database, perform BLAST searches to associate protein sequences to Metacyc reactions
```bash
diamond blastp --query query_protein.fasta --db  uniprot-seq-ids.db.dmnd --out query_protein_top5hit.blst \
               --outfmt 6 --evalue 0.001 --max-target-seqs 5 --sensitive

```
### step2 Summarize the reactions mapped to each query protein

```bash
python creat_RXN_dictionary.py query_protein_top5hit.blst protein_RXN_dictionary
```
link uniprot ID with Metacyc Reaction IDs, and creat dictionary with each query protein as the key and the RXNs for this 
protein as the values. Actually, here we did not use `uniprot-seq-ids.dat` to link uniprot ID with Metacyc Reaction IDs,
we just use information from `uniprot-seq-ids.fasta`.

### step3 associate Metacyc reactions to Metacyc pathways

According to `All_pathways_of_MetaCyc.txt`, map Reaction IDs to pathways.
