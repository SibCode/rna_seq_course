# rna_seq_course
RNA Sequencing Course Bioinformatics Master University of Bern

Folder structure is as follows:
scripts (contains all scripts)
- 01_run_fastqc.sh
- 02_get_references.sh
- ... etc. (for each step)
logs (contains all SLURM error and output logs)
- 01_run_fastqc
- 02_get_references
- ... etc. (for each script)
outputs (contains all output files of scripts)
- 01_run_fastqc
- 05_run_hisat2
- ... etc. (for each script that has outputs)
metadata (contains all metadata used by / for scripts)
references (contains reference genome for samples)

The last output folder is used as input for the next script (usually).

Due to large data the outputs, logs, and references folder is .gitignored.

Scripts were run in order.

Samplelist was created using the code from 04_get_samplelist.sh to create the file "samplelist.tsv" in the metadata folder.
