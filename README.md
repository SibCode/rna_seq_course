# rna_seq_course
RNA Sequencing Course Bioinformatics Master University of Bern

Folder structure is as follows:

- scripts (contains all scripts)
  - 01_run_fastqc.sh
  - 02_get_references.sh
  - ... etc. (for each step)
- logs (contains all SLURM error and output logs)
  - 01_run_fastqc
  - 02_get_references
  - ... etc. (for each script)
- outputs (contains all output files of scripts)
  - 01_run_fastqc
  - 05_run_hisat2
  - ... etc. (for each script that has outputs)
- metadata (contains all metadata used by / for scripts)
- references (contains reference genome for samples)

Due to data amount outputs, logs, and references folder is .gitignored.

Scripts were run in order unless specified. 10_Differential_Gene_Expression_Analysis.R contains the full analysis.

Samplelist was created using the code from 04_get_samplelist.sh to create the file "samplelist.tsv" in the metadata folder.
