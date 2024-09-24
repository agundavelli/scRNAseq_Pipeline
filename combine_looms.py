import os
import loompy

# Directories
loom_dir = "/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/final_results/final_looms"
combined_loom_dir = "/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/final_results/final_combined_looms"

# List of Loom file paths
looms = [
    "1D-GFP.loom", "1D-mCherry.loom", "1W-HFD.GFP.loom", 
    "1W-HFD.mCherry.loom", "2D-GFP.loom", "2D-mCherry.loom",
    "2W-HFD.GFP.loom", "2W-HFD.mCherry.loom", "3D-GFP.loom",
    "3D-mCherry.loom", "7D-GFP.loom", "7D-mCherry.loom",
    "RT-GFP.loom", "RT-mCherry.loom"
]

# List of specific loom file paths
GFP_looms = [
    "1D-GFP.loom", 
    "1W-HFD.GFP.loom", 
    "2D-GFP.loom",
    "2W-HFD.GFP.loom",
    "3D-GFP.loom",
    "7D-GFP.loom", 
    "RT-GFP.loom"
]

mCherry_looms = [
    "1D-mCherry.loom", 
    "1W-HFD.mCherry.loom", 
    "2D-mCherry.loom",
    "2W-HFD.mCherry.loom",
    "3D-mCherry.loom",
    "7D-mCherry.loom",
    "RT-mCherry.loom"
]

HFD_looms = [
    "1W-HFD.GFP.loom",
    "2W-HFD.GFP.loom",
    "1W-HFD.mCherry.loom",
    "2W-HFD.mCherry.loom",
    "RT-GFP.loom",
    "RT-mCherry.loom"
]

cold_looms = [
    "1D-GFP.loom",
    "2D-GFP.loom",
    "3D-GFP.loom",
    "7D-GFP.loom", 
    "RT-GFP.loom",
    "1D-mCherry.loom",  
    "2D-mCherry.loom",
    "3D-mCherry.loom",
    "7D-mCherry.loom",
    "RT-mCherry.loom"
]

def combine_loom_files(loom_dir, combined_looms_dir, loom_files, output_filename, key):
    loom_file_paths = [os.path.join(loom_dir, file) for file in loom_files]
    
    # Combine loom files
    loompy.combine(loom_file_paths, output_file=os.path.join(combined_looms_dir, output_filename), key=key)

# List of all combinations
loom_combinations = [
    (looms, "all_combined.loom"),
    (GFP_looms, "GFP_combined.loom"),
    (mCherry_looms, "mCherry_combined.loom"),
    (HFD_looms, "HFD_combined.loom"),
    (cold_looms, "cold_combined.loom")
]

# Run the combine function for each combination
for loom_files, output_filename in loom_combinations:
    combine_loom_files(loom_dir, combined_loom_dir, loom_files, output_filename, key="Accession")
