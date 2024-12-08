#!/usr/bin/env bash
#SBATCH --job-name=feature_count
#SBATCH --output=/data/users/sbertschinger/rna_seq_course/logs/09_run_feature_count/output_%j.out
#SBATCH --error=/data/users/sbertschinger/rna_seq_course/logs/09_run_feature_count/error_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=simon.bertschinger@students.unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8

# Note: Slurm parameter array is set to 1-16 for the 16 samples
WORKDIR="/data/users/sbertschinger/rna_seq_course"
INDIR="$WORKDIR/outputs/07_run_samtools_sort"
OUTDIR="$WORKDIR/outputs/09_run_feature_count"

# Define input BAM file and output BAM file
INPUT_FILES="${INDIR}/*.bam"
OUTPUT_FILE="${OUTDIR}/feature_counts.txt"

# Define reference genome used
REFERENCE_FILE="$WORKDIR/references/Mus_musculus.GRCm39.113.gtf"

# Log output file for reference in case of sample specific error
echo "Run feature count for $INPUT_FILES..."

# Run the apptainer feature count with a minimum mapping quality score of 10 (-Q 10)
apptainer exec --bind /data/ /containers/apptainer/subread_2.0.1--hed695b0_0.sif\
    featureCounts ${INPUT_FILES} -p -s 2 -T $SLURM_CPUS_PER_TASK -t exon -g gene_id\
    -a ${REFERENCE_FILE} -o ${OUTPUT_FILE} -Q 10

# Log succesful feature counts run for file with output file path
echo "Ran feature count for '$INPUT_FILES' and saved output file to: $OUTPUT_FILE"
