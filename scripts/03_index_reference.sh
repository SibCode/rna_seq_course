#!/usr/bin/env bash
#SBATCH --job-name=index_reference
#SBATCH --output=/data/users/sbertschinger/rna_seq_course/logs/03_index_reference/output_%j.out
#SBATCH --error=/data/users/sbertschinger/rna_seq_course/logs/03_index_reference/error_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=simon.bertschinger@students.unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=pibu_el8

PATH_TO_OUTPUT=/data/users/sbertschinger/rna_seq_course/references

apptainer exec --bind /data/ \
    /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
    hisat2-build -p 4 ${PATH_TO_OUTPUT}/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    ${PATH_TO_OUTPUT}/Mus_musculus.GRCm39.dna.primary_assembly
