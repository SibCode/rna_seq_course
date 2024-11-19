#!/usr/bin/env bash
#SBATCH --job-name=hisat2_array
#SBATCH --output=/data/users/sbertschinger/rna_seq_course/logs/05_run_hisat2/output_%j.out
#SBATCH --error=/data/users/sbertschinger/rna_seq_course/logs/05_run_hisat2/error_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=simon.bertschinger@students.unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=1-16

# Note: Slurm parameter array is set to 1-16 for the 16 samples
WORKDIR="/data/users/sbertschinger/rna_seq_course"
OUTDIR="$WORKDIR/outputs/05_run_hisat2"
SAMPLELIST="$WORKDIR/metadata/samplelist.tsv"
REFERENCE_PATH="$WORKDIR/references"

# Get the sample name out of the samplelist, as well as the two reads
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# Define output file
OUTPUT_FILE="${OUTDIR}/${SAMPLE}.sam"

# Log above in the output file for reference in case of an error
echo "Run hisat2 for $SAMPLE with $READ1 and $READ2 ..."

# Run the apptainer hisat2 using rna-strandness RF as samples are paired-end and strand specific
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif\
    hisat2 -x ${REFERENCE_PATH}/Mus_musculus.GRCm39.dna.primary_assembly\
    -1 ${READ1} -2 ${READ2} -S ${OUTPUT_FILE} -p $SLURM_CPUS_PER_TASK --rna-strandness RF
