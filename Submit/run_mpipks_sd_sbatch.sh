#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00

scratch=/scratch/ivanchen/kalyakulina/mtDNA/$1
code_base=/home/ivanchen/kalyakulina/mtDNA/Source/mtDNA/cluster
mkdir -p $scratch
mkdir -p $1
cd $scratch
cp $1/config.txt .
cp $1/config_mt_genes.txt .
cp $1/config_nuc_genes.txt .

cat config.txt

srun python $code_base/random_forest.py

cp -r $scratch/* $1 # Better use a subdirectory of $HOME .
rm -r $scratch/*
