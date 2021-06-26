# Metagenomics workflow from raw sequencing data to metagenome-assembled genomes (MAGs)

This workflow is generated to process the analysis for metagenomics sequencing samples in TEDDY
study (dbGaP Study Accession: phs001443.v1.p1).

**Data processing for metagenomics assembly**
  * a. Tools
  * b. Download and PreProcess the Raw Data
  * c. Assembly
  * d. Quantification of scaffolds
  * e. ORF prediction and functional annotation with KEGG and MetaCyc
  * f. Genome binning and Quantification of MAGs


## Tools
1. NCBI SRA Toolkit to download the sequencing file with sra format, and then convert it to fastq format
2. [**BBtools**]:https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/ for raw data preprocessing 
3. Assembly tool **metaspades** https://github.com/ablab/spades
4. **pullseq** to filter sequences by a minimum length or maximum length https://github.com/bcthomas/pullseq 
5. reads mapping and calculation
**bowtie** for read mappping https://github.com/BenLangmead/bowtie2
**shrinksam** to remove sequences that failed to map to the reference https://github.com/bcthomas/shrinksam
**add_read_count.rb** https://github.com/bcthomas/misc_scripts/tree/master/add_read_count
5. **Prodigal** for ORF prediction https://github.com/hyattpd/Prodigal
6. **KofamScan** KEGG ko term annotation  https://github.com/takaram/kofam_scan
7. **Diamond**  alignment https://github.com/bbuchfink/diamond

## Download and PreProcess the Raw Data

**1. download data from dbGaP**

download data and covert to fastq format using NCBI SRA Toolkit
```bash
$ prefetch --option-file sra_list.txt

# an example for `sra_list.txt`:
SRR7560855
SRR7560856
SRR7560857
```
**2. convert sra files to fastq format**

In the following assembly process, we will co-assemly the samples from the same subject.
Therefore, we merged the fastq files from the same subject.
```bash
$ for i in *.sra
  do
    fastq-dump --split-files --gzip $i
  done

$ cat *_1.fastq.gz >SUBJECTID_1.fastq.gz
$ cat *_2.fastq.gz >SUBJECTID_2.fastq.gz
```
**3. PreProcess the Raw Data**

The metagenomics data from this dbGaP Study (phs001443.v1.p1) is well proprocessed with quality control and removing human genome sequences.
If your data is the raw sequencing data without preprocess, BBtools could be used to do this step.
## Assembly

**1. metaspades**
```bash
$ metaspades.py -1 /work/TEDDY/ncbi/dbGaP-21828/sra/SUBJECTID/SUBJECTID_1_.fastq.gz \
                -2 /work/TEDDY/ncbi/dbGaP-21828/sra/SUBJECTID/SUBJECTID_2_.fastq.gz \
                -o SUBJECTID \
                -t 20  
```
The default memorial is 250G, if you need a larger one, please add a parameter setting like -m 500 to 500G.

**2. N50, N90 Calculation (BBtools)**
```bash
$ stats.sh in=assembly_scaffolds.fasta out=assembly_data_stats
```
## Abundance Calculation for each scaffold
**1. fillter scaffolds**
```bash
#Select scaffolds that sequence length > 1000; 
$ pullseq -i SUBJECTID_scaffolds.fasta -m 1000 > SUBJECTID_scaffolds_min1000.fasta
```
**2. map reads into scaffolds with bowtie2**
```bash
#Create bowtie2 index file
$ bowtie2-build  SUBJECTID_scaffolds_min1000.fasta  SUBJECTID_min1000

# Sequence alignment for each SRR runs; the reference is the co-assemblied scaffolds of the subject where SRR runs from.

$ bowtie2 -x SUBJECTID_min1000 -1 SUBJECTID_eachSRR_1.fastq.gz \
                               -2 SUBJECTID_eachSRR_2.fastq.gz \
                               -S SUBJECTID_eachSRR_alignment.sam \
                               -p 19

# 'SUBJECTID_eachSRR_alignment.sam' is also the file used for genome binning
# You can obtain the number of aligned reads in the output file 
# and the number of total read can be obtained in the log file of assembly
# Mapping rate = the number of aligned reads/ the number of total readsRPKM Calculation (shrinksam)
```
**3. Count Mapped Reads**
```bash
# we use `shrinksam` to remove unmapped reads from bowtie2-created SAM files,
# it will generate a SAM file with only the reads that map to an assembled scaffold.

$ shrinksam -i SUBJECTID_eachSRR_alignment.sam > SUBJECTID_eachSRR_mapped.sam

# Then we count the reads numbers mapped to each scaffold through `add_read_count.rb`
$ add_read_count.rb SUBJECTID_eachSRR_mapped.sam  SUBJECTID_scaffolds_min1000.fasta > SUBJECTID_eachSRR.reads.counted

# we just filter the lines containing scaffold name in the output fasta files

$ grep -e ">" SUBJECTID_eachSRR.reads.counted > mapped.counted.result
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
$ Python parse_scaffold_rpkm.py mapped.counted.result total_reads_file srrID

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

## Protein prediction and read counts for each protein (ORF) 
1. Protein prediction through `Prodigal`

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
`We divide the scaffold read counts into each protein according to their CDS length`. 


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
For more information: https://github.com/thepanlab/Seq2MAG/blob/master/annotation/metacyc/metacyc_annotation.md

## Genome binning
**1.sam to sorted bam for binning**
```bash
$ samtools view -bS SUBJECTID_eachSRR_alignment.sam -@ 19 > SUBJECTID_eachSRR_alignment.bam
$ samtools sort SUBJECTID_eachSRR_alignment.bam  -@ 19 -o SUBJECTID_eachSRR_sorted > SUBJECTID_eachSRR_sorted.bam

# SUBJECTID_eachSRR_alignment.sam is from Bowtie2 output result in read mapping process.
```
**2. metabat binning**
```bash
$ jgi_summarize_bam_contig_depths --outputDepth output_depth.txt *_sorted.bam
# *_sorted.bam are the bam files from the same subjects.

$ metabat -i sequence_min1000.fasta -o output_bin -a output_depth.txt -m 2000
```
**3. qualtify for the binned genomes**
```bash
$ checkm lineage_wf -f CheckM.txt -t 10 -x fa --pplacer_threads 1 input_fold_includes_Bins checkm_wf_out  
```




