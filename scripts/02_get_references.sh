#!/usr/bin/env bash
#SBATCH --job-name=get_references
#SBATCH --output=/data/users/sbertschinger/rna_seq_course/logs/02_get_references/output_%j.out
#SBATCH --error=/data/users/sbertschinger/rna_seq_course/logs/02_get_references/error_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=simon.bertschinger@students.unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8

PATH_TO_OUTPUT=/data/users/sbertschinger/rna_seq_course/references

wget https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz $PATH_TO_OUTPUT
wget https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz $PATH_TO_OUTPUT

sum ${PATH_TO_OUTPUT}/*

# Unzip the files
gunzip ${PATH_TO_OUTPUT}/*.gz