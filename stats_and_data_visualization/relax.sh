#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name="relax"
#SBATCH --mem 100G
#SBATCH --output=relax_02.log
#SBATCH --exclude=node117,node118

module purge
module load anaconda/colsa
#conda activate hyphy-2.5.26
conda activate hyphy-2.5.59

for x in /mnt/lz01/plachetzki/hkp1001/chapter2/concatenated_files/all_labeled/*.fa; \
do hyphy relax --alignment ${x} --test Foreground;
done

