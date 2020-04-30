# Metagenomics workflow from raw sequencing data to metagenome-assembled genomes (MAGs)

1. Preprocessing of raw sequencing data
2. Metagenome assembly
3. Quantification of scaffolds
4. Functional annotation of proteins with UniFam and KEGG
5. Binning and quantification of MAGs
6. Taxanomy analysis of MAGs
7. MetaCyc pathway reconstruction of MAGs

7.1  using 'uniprot-seq-ids.fasta' as blast database, perform BLAST searches to associate protein sequences to Metacyc reactions

**A common early step in performing pathway analysis of genomes and metagenomes is to associate protein sequences to MetaCyc reactions.** The Pathway Tools software infers such associations by using EC numbers, enzyme names, and Gene Ontology terms within protein annotations. Such annotations might be inferred using a variety of sequence-analysis methods.

To aid researchers in associating sequences to MetaCyc reactions, MetaCyc enzymes that have a link to UniProt contain protein sequence information. It is possible to perform BLAST searches against MetaCyc proteins with sequence information using the "BLAST Search" command under the Search menu. In addition, each release of MetaCyc includes **a file that associates MetaCyc reaction IDs with the UniProt identifiers of enzymes known to catalyze those reactions. Note that not all MetaCyc reactions have EC numbers (because not all enzymes have yet been assigned EC numbers), therefore EC numbers are not a comprehensive mechanism for associating sequences to reactions.** The file is called **uniprot‑seq‑ids.dat** and is included in the MetaCyc data file distribution.
if you want to know how to get the uniprot-seq-ids.dat, please see details in the following links:
https://metacyc.org/MetaCycUserGuide.shtml;
http://bioinformatics.ai.sri.com/ptools/flatfile-format.html

Using DIAMOND to perform blast searches:
```
diamond blastp --query SUBJECTID_protein.fasta --db  uniprot-seq-ids.db.dmnd --out SUBJECTID_protein_top5hit.blst --outfmt 6 --evalue 0.001 --max-target-seqs 5 --sensitive

python creat_RXN_dictionary.py SUBJECTID_protein_top5hit.blst SUBJECTID_RXN_dic
```

7.2 associate Metacyc reactions to Metacyc pathways
login to Metacyc database, and then creat your own specical smartable of 'All pathways of MetaCyc' (you can add columns in this samrtable to get more information of each pathway and download this smartable). 
https://metacyc.org/;
https://metacyc.org/PToolsWebsiteHowto.shtml#node_sec_6
