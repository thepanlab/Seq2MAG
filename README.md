# Metagenomics workflow from raw sequencing data to metagenome-assembled genomes (MAGs)

This pipline is generated to process the metagenomics analysis for samples in TEDDY
study.

**Data processing for metagenomics assembly**
  * a. Tools
  * b. PreProcess the Raw Data
  * c. Assembly
  * d. Abundance calculation
  * e. ORF prediction and functional annotation
  * f. Genome binning
## Tools
1. [**BBtools**]:https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/ for raw data preprocessing 
2. Assembly tools **metaspades** https://github.com/ablab/spades or **megahit**  https://github.com/voutcn/megahit
3. **pullseq** to filter sequences by a minimum length or maximum length https://github.com/bcthomas/pullseq 
4. reads mapping and calculation
**bowtie** for read mappping https://github.com/BenLangmead/bowtie2
**shrinksam** to remove sequences that failed to map to the reference https://github.com/bcthomas/shrinksam
**add_read_count.rb** https://github.com/bcthomas/misc_scripts/tree/master/add_read_count
5. **Prodigal** for ORF prediction https://github.com/hyattpd/Prodigal
6. **KofamScan** KEGG ko term annotation  https://github.com/takaram/kofam_scan
7. **Diamond**  alignment https://github.com/bbuchfink/diamond
8. **Deseq** Statistical Analysis ttps://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf

## PreProcess the Raw Data (BBtools)

**1. Trim adaptors**

Besides trim adaptors, optionally, reads with Ns can be disregarded by adding “maxns=0” ;
and reads with really low average quality can be disregarded with “maq=8”
```bash
$ bbduk.sh in=input.fastq out=output.fastq \
                          ref=${BBtools_path}/resources/adapters.fa \
                          ktrim=r k=23 mink=11 hdist=1 tpe tbo
```
**2. Remove synthetic artifacts and spike-ins by kmer-matching**

This will remove all reads that have a 31-mer match to PhiX, allowing one mismatch.
The “outm” stream will catch reads that matched a reference kmers;
“stats” will produce a report of which contaminant sequences were seen, and how many reads had them.
```bash
$ bbduk.sh in=trimed.fastq out=output_u.fastq outm=output_m.fastq \
           ref=${BBTools_path}/resources/sequencing_artifacts.fa.gz \
           ref=${BBTools_path}/resources/phix174_ill.ref.fa.gz \
           k=31 hdist=1 stats=stats.txt
```
**3. Error Correction**

Low-depth reads can be discarded here with the”tossdepth”, or “tossuncorrectable” flags
For very large datasets, “prefilter=1” or “prefilter=2” can be added to conserve memory
```bash
$ tadpole.sh -Xmx32g in= output_u.fastq out=output_tecc.fastq filtermemory=7g ordered prefilter=1
```
## Assembly (Metaspades or Megahit)

Metaspades require more memorial and longer time, but can get a higher quality assembly and scaffold directly
Megahit require less memorial and less time, but with relatively lower quality and can only assembly the reads into contig.

**1. metaspades**
```bash
$ metaspades --12  output_tecc.fastq -o  assembly_data.fasta
#--12 means the sequencing data is in interleaved form
```

The default memorial is 250G, if you need a larger one, you can set like -m 500

**2. Megahit**
```bash
$ megahit --12  output_tecc.fastq -o assembly_data_megahit.fasta
```

**3. N50, N90 Calculation (BBtools)**
```bash
$ stats.sh in=assembly_data.fasta out=assembly_data_stats
```
## Abundance Calculation for each scaffold
**1. fillter scaffolds**
```bash
#Select scaffolds that sequence length > 1000; 
$ pullseq -i assembly_data.fasta -m 1000 > sequence_min1000.fasta
```
**2. map reads into scaffolds with bowtie2**
```bash
#Create bowtie2 index file
$ bowtie2 -build  sequence_min1000.fasta sequence_min1000

#Deinterleave paired reads
$ reformat.sh in=output_tecc.fastq out1=read1.fastq out2=read2.fastq

#Sequence alignment
$ bowtie2 -x sequence_min1000 -1 read1.fastq -2 read2.fastq S alignment.sam -p 19

# 'alignment.sam' is also the file used for genome binning

# You can obtain the number of aligned reads in the output file 
# and the number of total read can be obtained in the log file of assembly
# Mapping rate = the number of aligned reads/ the number of total readsRPKM Calculation (shrinksam)
```
**3. Count Mapped Reads**
```bash
# we use `shrinksam` to remove unmapped reads from bowtie2-created SAM files,
# it will generate a SAM file with only the reads that map to an assembled scaffold.

$ shrinksam -i alignment.sam > mapped.sam

# Then we count the reads numbers mapped to each scaffold through `add_read_count.rb`
$ add_read_count.rb mapped.sam sequence_min1000.fastq > mapped.fasta.counted

# we just filter the lines containing scaffold name in the output fasta files

$ grep -e ">" mapped.fasta.counted > mapped.counted.result
$ sed -i "s/>//g" mapped.counted.result
$ sed -i "s/read_count_//g" mapped.counted.result

# an example of `mapped.count.result`
NODE_8_length_295416_cov_86.298909  read_length_150    1805
NODE_9_length_294755_cov_77.048012  read_length_150     267
NODE_10_length_287310_cov_147.7697  read_length_150     2469
NODE_11_length_283217_cov_82.5382   read_length_150     423
NODE_12_length_273911_cov_216.2412  read_length_150     16

```
**4. RPKM Calculation for each assembled scaffold**
```bash
$ Python parse_rpkm.py mapped.counted.result total_reads_file srrID

# total_reads_file is a file that includes total reads number in each sample, an example:
10004483	SRR7768473
10005994	SRR7768567
10007320	SRR7768664

# srrID is just the srrID sample ID number (such as SRR7768473)

# to calculate the rpkm, scaffold length is all needed. 
# MetaSpades predicted scaffold name is like 'NODE_8_length_295416_cov_86.298909'.
# the number 295416 is the length of this scaffold.
```
**Important notes: Here we only got RPKM for each scaffold; 
For the proteins in the same scaffod, they should have the same RPKM with this scaffold.**

## ORF Prediction and read counts for each ORF 
1. ORF prediction through `Prodigal`

```bash
$ prodigal -i sequence_min1000.fasta -p meta -a trans_protein.fasta -f gff -o predicted.gff

# two output files: 
trans_protein.fasta (protein translations file); 
predicted.gff (inforamtion for each CDS, we will use the CDS length) 
```
**2.1 Map scaffold RPKM to each protein that located in this scaffold.**

The principle is that for the proteins in the same scaffod, they should have the same RPKM with this scaffold. 
This is analysis for metagenome, which is not same with metatranscriptomics.

**2.2 read counts for each protein**

Some statistic tools, such as DEGseq2, they use the read counts (not the normalized RPKM) as the inputs
here we add a method to get protein read counts directly from the scaffold read counts in the above.
`We could just divide the scaffold read counts into each protein according to their CDS length`. 


## Functional Annotation of the predicted ORF

**1.	KEGG annotation through `KofamScan`**
```bash
$ kofam_scan/exec_annotation -o Coassembly_KO.txt trans_protein.fasta --tmp-dir tmp_KO --cpu 10
```
map the annotated KO term into pathways through `KEGG Mapper`: https://www.genome.jp/kegg/tool/map_pathway.html

**2. Metacyc reaction annotation by aligning database from Metacyc pathway database;
`diamond` was used to do sequence alignment**
```bash
$ diamond blastp --query trans_protein.fasta \
                 --db uniprot-seq-ids.fasta \
                 --out Metacyc_protein_top5hit.blst \
                 --outfmt 6  --evalue 0.001 --max-target-seqs 5 --sensitive
                 
$ python creat_RXN_dictionary.py Metacyc_protein_top5hit.blst Metacyc_protein_RXN_key_sen
```
For more information: https://github.com/thepanlab/MetagenomicsTree/blob/master/metacyc/metacyc_annotation.md

## Genome binning
**1.sam to sorted bam for binning**
```bash
$ samtools view -bS alignment.sam -@ 19 > alignment.bam
$ samtools sort alignment.bam  -@ 19 -o alignment_sorted > alignment_sorted.bam

# alignment.sam is from Bowtie2 output result in read mapping process.
```
**2. metabat binning**
```bash
$ jgi_summarize_bam_contig_depths --outputDepth output_depth.txt *_sorted.bam
$ metabat -i sequence_min1000.fasta -o output_bin -a output_depth.txt -m 2000
```
**3. qualtify and taxonomy annotation for the binned genomes**
```bash
$ checkm lineage_wf -f CheckM.txt -t 10 -x fa --pplacer_threads 1 input_fold_includes_Bins checkm_wf_out  
```






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
