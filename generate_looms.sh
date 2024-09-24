#!/bin/bash
#$ -cwd
#$ -P farmer
#$ -pe omp 28
#$ -l mem_per_core=13G
#$ -l h_rt=3:00:00  # Adjust the walltime as needed based on your job's expected duration
#$ -N generate_looms_RT-mCherry

# Load necessary modules
# Make sure you conda activate loom_env
#module load samtools/1.12
#module load velocyto

# Directories and file paths
bam_dir="/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/final_results/bam_files"
barcodes_dir="/project/farmer/samples"
results_dir="/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/final_results/final_looms"
gtf_file="/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/Mus_musculus.GRCm39.112.gtf"

condition="RT-mCherry"
bam_file="${bam_dir}/${condition}.ScaleRNA.bam"
barcode_file="${barcodes_dir}/${condition}.ScaleRNA.filtered.matrix/barcodes.tsv"

# Log file path (adjust as per your preference)
log_file="/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/logs/sort_and_velocyto_RT-mCherry.log"

echo "Condition: ${condition}"
echo "BAM file: ${bam_file}"
echo "Barcode file: ${barcode_file}"

# Redirect stdout and stderr to log file
exec >"${log_file}" 2>&1

# Run Velocyto
echo "Running Velocyto for ${condition}..."
velocyto run -b "${barcode_file}" -o "${results_dir}" "${bam_file}" "${gtf_file}"

echo "Velocyto processing for ${condition} completed."
