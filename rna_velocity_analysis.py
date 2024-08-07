import os
import loompy
import anndata as ad

import warnings
warnings.simplefilter("ignore", category=UserWarning)

import scanpy as sc
import scvelo as scv
import cellrank as cr
import matplotlib.pyplot as plt

# Adjust Scanpy and scVelo settings if needed
sc.settings.set_figure_params(frameon=False, dpi=100)
cr.settings.verbosity = 2

# Directories
loom_dir = "/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/results/loom_files"
#loom_10x_dir = "/projectnb/farmer/akhila/10X-Rabhi/"
results_dir = "/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/results" 
#results_10x_dir = "/projectnb/farmer/akhila/10X-Rabhi/results"

# Path to the combined Loom file
combined_loom_path = "/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/results/loom_files/cellsorted_1D-GFP.loom"
#loom_path = "/projectnb/farmer/akhila/10X-Rabhi/FS_veh.loom"

# Load the Loom file
adata = sc.read_loom(combined_loom_path)

#############################################################################################################################
# Step 1: Compute and visualize the fraction of spliced and unspliced reads
#############################################################################################################################
# Check layers in the AnnData object
print("Available layers in the AnnData object:")
print(adata.layers.keys())

print("Computing proportions.")

# Calculate proportions of unsliced and spliced reads
scv.pl.proportions(adata)

# Save the figure 
plt.savefig(os.path.join(results_dir, "proportions.png"))

# Show the figure
plt.show()
plt.close()

print("Proportions plot saved.")

#############################################################################################################################
# Step 2: Compute RNA velocity
#############################################################################################################################

# Read the Loom file
#adata = scv.read(combined_loom_path, cache=True)

# Compute PCA
print("Computing PCA.")
sc.pp.pca(adata)

# Compute neighbors
print("Computing neighbors.")
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)

# Compute UMAP embedding
print("Computing UMAP embedding.")
sc.tl.umap(adata)

# Compute moments
print("Computing moments.")
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

# Compute RNA velocity
print("Computing RNA velocity.")
scv.tl.velocity(adata)

# Compute RNA velocity graph
scv.tl.velocity_graph(adata)

# Run scVelo Analysis
print("Running scVelo analysis (recover dynamics).")
scv.tl.recover_dynamics(adata, n_jobs=8) # Recover dynamics
print("Computing RNA velocity (dynamical mode).") # Compute RNA velocity (dynamical mode)
scv.tl.velocity(adata, mode="dynamical") # Recompute velocity graph
print("Recomputing velocity graph.") 
scv.tl.velocity_graph(adata)

# Plot RNA velocity
print("Plotting RNA velocity.")

# Check and list embeddings in the AnnData object
print("Available embeddings in the AnnData object:")
print(adata.obsm.keys())

# Check available layers
print(adata.layers)

'''
# Load your AnnData object
adata_file = os.path.join(results_dir, "loom_with_velocity.h5ad")
adata = sc.read(adata_file)
'''

# Plot RNA velocity with streamlines and a legend on the right
fig, ax = plt.subplots(figsize=(10, 8))  # Adjust the figure size as needed
scv.pl.velocity_embedding_stream(adata, basis='umap', color='CellType', ax=ax, legend_loc='right', legend_fontsize=12)
#legend_loc='right', legend_fontsize=12

# Adjust the legend position
legend = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()

# Save the figure
plt.savefig(os.path.join(results_dir, 'velocity_streamline_plot.png'), bbox_inches='tight')

# Show the figure
plt.show()
plt.close()

print("Streamline velocity plot saved.")

# Save the AnnData object
adata_file = os.path.join(results_dir, "adata_with_velocity.h5ad")
adata.write(adata_file)
print(f"AnnData object saved to {adata_file}")
