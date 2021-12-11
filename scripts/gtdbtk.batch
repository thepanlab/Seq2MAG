#!/bin/bash
#SBATCH --partition=omicsbio
#SBATCH --ntasks=1
#SBATCH --nodelist=c660
#SBATCH --time=200:00:00
#SBATCH --job-name=GTDB
#SBATCH --mem=180G
#SBATCH --workdir=/work/TEDDY/binning/GTDB/
#SBATCH --output=GTDB_%J_stdout.txt
#SBATCH --error=GTDB_%J_stderr.txt
#SBATCH --mail-user=lizhang12@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/work/TEDDY/binning/GTDB/

module load Python/3.6.3-intel-2016a
module load HMMER/3.2.1-foss-2018b
module load pplacer/1.1.alpha19


#wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
#tar xvzf gtdbtk_r89_data.tar.gz

gtdbtk classify_wf --batch TEDDY_MAGs --out_dir TEDDY_GTDB --cpus 40 