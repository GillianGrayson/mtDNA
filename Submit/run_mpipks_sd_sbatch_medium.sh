#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=5000
#SBATCH --time=48:00:00
#SBATCH --partition=medium

scratch=/scratch/denysov/yusipov/mtDNA/$1
code_base=/home/denysov/yusipov/mtDNA/Source/mtDNA/cluster
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
