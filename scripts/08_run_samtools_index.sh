#!/usr/bin/env bash
#SBATCH --job-name=samtools_index_array
#SBATCH --output=/data/users/sbertschinger/rna_seq_course/logs/08_run_samtools_index/output_%j.out
#SBATCH --error=/data/users/sbertschinger/rna_seq_course/logs/08_run_samtools_index/error_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=simon.bertschinger@students.unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=1-16

# Note: Slurm parameter array is set to 1-16 for the 16 samples
WORKDIR="/data/users/sbertschinger/rna_seq_course"
INDIR="$WORKDIR/outputs/07_run_samtools_sort"
OUTDIR="$WORKDIR/outputs/08_run_samtools_index"
SAMPLELIST="$WORKDIR/metadata/samplelist.tsv"

# Get the sample name out of the samplelist
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`

# Define input SAM file and output BAM file
INPUT_FILE="${INDIR}/${SAMPLE}.bam"
OUTPUT_FILE="${OUTDIR}/${SAMPLE}.bam"

# Log output file for reference in case of sample specific error
echo "Run samtools index for $SAMPLE ..."

# Run the apptainer samtools index task and output in specific folder
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif\
    samtools index ${INPUT_FILE} -o ${OUTPUT_FILE}
