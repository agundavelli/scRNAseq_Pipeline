#!/bin/bash
#$ -cwd
#$ -P farmer
#$ -pe omp 28
#$ -l mem_per_core=18G
#$ -l h_rt=4:00:00  # Adjust the walltime as needed based on your job's expected duration
#$ -N moscot

# use pe omp 28 and mem_per_core-18G
# Load required modules if necessary
#module load python
#module load miniconda

# Activate your Python environment if not already activated
#conda activate rabhi

# Run your Python script
#python CellRank2_Analysis.py
#python generate_heatmap.py
#python combine_looms.py
python rna_velocity_analysis.py
#python add_metadata_to_looms.py