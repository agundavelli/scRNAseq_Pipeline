import sys
import os

# Activate moscot environment 
# Import necessary modules
from moscot.problems.time import TemporalProblem
import cellrank as cr
import scanpy as sc
from cellrank.kernels import RealTimeKernel
import matplotlib.pyplot as plt

# Set Scanpy and CellRank settings
sc.settings.set_figure_params(frameon=False, dpi=100)
cr.settings.verbosity = 2

import warnings
warnings.simplefilter("ignore", category=UserWarning)

# Results directory path
results_dir = "/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/final_results"

# Set Scanpy's figure directory to the results directory
sc.settings.figdir = results_dir


# Read in AnnData object
adata_path = os.path.join(results_dir, "rna_velocity", "all_combined_adata_with_velocity.h5ad")
print(f"Reading AnnData object from: {adata_path}")
adata = sc.read_h5ad(adata_path)
print(f"AnnData object loaded with {adata.n_obs} observations and {adata.n_vars} variables")

# Define a mapping for 'day' column to numerical values
day_mapping = {
    "RT": 0,      # Room Temperature
    "1D": 1,      # 1 Day
    "2D": 2,      # 2 Days
    "3D": 3,      # 3 Days
    "7D": 7,      # 7 Days
    "1W": 7,      # 1 Week (equivalent to 7 days)
    "2W": 14      # 2 Weeks (equivalent to 14 days)
}

# Convert the 'day' column to categorical and numerical types
print("Converting 'day' column to categorical and numerical types")
adata.obs["day"] = adata.obs["day"].astype("category")
adata.obs["day_numerical"] = adata.obs["day"].map(day_mapping).astype(float)

# Filter the data for "cold" and "RT" conditions only
conditions_of_interest = ["RT", "1D", "2D", "3D", "7D"]  # Adjust based on your specific cold conditions
#conditions_of_interest = ["RT", "1W", "2W"]  # Adjust based on your specific HFD conditions
print(f"Filtering data for conditions: {conditions_of_interest}")
adata_filtered = adata[adata.obs["day"].isin(conditions_of_interest), :]
print(f"Filtered data now has {adata_filtered.n_obs} observations")

# Optionally, subsample the filtered data to speed up the analysis
print("Subsampling data to 25% of its original size")
sc.pp.subsample(adata_filtered, fraction=0.25, random_state=0)
print(f"Subsampled data now has {adata_filtered.n_obs} observations")

# Define file paths for saving plots
umap_plot_filename = "day_celltype_cold_RT.png"
# Visualize the filtered data to get an understanding of its distribution
print("Visualizing data with UMAP embeddings")
sc.pl.embedding(
    adata_filtered,
    basis="umap",  # or use another embedding if you prefer
    color=["day_numerical", "CellType"],  # Update with appropriate metadata columns
    color_map="gnuplot",
    save=umap_plot_filename  # Save the plot as a PNG file in the specified directory
)

# Compute PCA and the k-NN graph needed for the Moscot analysis
print("Computing PCA")
sc.pp.pca(adata_filtered)
print("Computing k-NN graph")
sc.pp.neighbors(adata_filtered, random_state=0)

# Initialize the TemporalProblem with your filtered AnnData object
print("Initializing TemporalProblem")
tp = TemporalProblem(adata_filtered)

# Adjust the marginals for cellular growth and death using appropriate gene sets
print("Scoring genes for marginals")
tp = tp.score_genes_for_marginals(
    gene_set_proliferation="mouse",
    gene_set_apoptosis="mouse"
)

proliferation_apoptosis_plot_filename = "proliferation_apoptosis_cold_RT.png"
# Visualize proliferation and apoptosis scores
print("Visualizing proliferation and apoptosis scores")
sc.pl.embedding(
    adata_filtered, basis="umap",  # or "force_directed" if you have that embedding
    color=["CellType", "proliferation", "apoptosis"],  # Update with your annotations
    save=proliferation_apoptosis_plot_filename,  # Save the plot as a PNG file in the specified directory
    size=30,  # Adjust point size
    alpha=0.6  # Adjust point transparency
)

tp = tp.prepare(time_key="day_numerical")
tp = tp.solve(epsilon=1e-3, tau_a=0.95, scale_cost="mean")

# Set up the RealTime kernel
tmk = RealTimeKernel.from_moscot(tp)
tmk.compute_transition_matrix(self_transitions="all", conn_weight=0.2, threshold="auto")

# Save the random walks plot
random_walks_plot_filename = os.path.join(results_dir, "cold_random_walks_plot.png")
print("Saving random walks plot")
tmk.plot_random_walks(
    max_iter=250,
    start_ixs={"day_numerical": 0.0},
    basis="force_directed",
    seed=0,
    dpi=100,
    size=30,
    save=random_walks_plot_filename  # Save the plot directly
)

# Save the flow plot between clusters (cell types)
flow_plot_filename = os.path.join(results_dir, "cold_flow_plot.png")
print("Saving flow plot")
ax = tmk.plot_single_flow(
    cluster_key="CellType",
    time_key="day_numerical",
    cluster="Endothelial Cells",
    min_flow=0.2,
    xticks_step_size=4,
    show=False,
    clusters=["M2 Macrophages", "M2 Tissue-Resident Macrophages", "M2/Tissue-Resident Macrophages", 
              "Subset Macrophages", "T Cells", "Immune Cells", "Neurons", 
              "Proliferating Cells", "Endothelial Cells", "Vascular Cells"],
)
_ = ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
plt.savefig(flow_plot_filename, dpi=300)
plt.close()  # Close the plot to free up memory

# Save the AnnData object
#new_adata_filename = os.path.join(results_dir, 'HFD_adata_with_transitions.h5ad')
#adata.write(new_adata_filename)
#print(f"AnnData object saved to: {new_adata_filename}")

# Delete adata_filtered to free up memory
del adata_filtered

print("Script completed successfully")