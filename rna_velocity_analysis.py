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
results_dir = "/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/final_results/rna_velocity" 
combined_loom_path = "/projectnb/farmer/akhila/Rabhi_Lab/scSPLiT-seq/final_results/copy_final_combined_looms/GFP_combined.loom"
# all_combined.loom, cold_combined.loom, GFP_combined.loom, HFD_combined.loom, mCherry_combined.loom

# Load the Loom file
adata = sc.read_loom(combined_loom_path)

# Save after loading the Loom file
#adata_file = os.path.join(result_dir, "loaded_adata.h5ad")
#adata.write(adata_file)
#print(f"AnnData object saved after loading Loom file to {adata_file}")

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
plt.savefig(os.path.join(results_dir, "GFP_combined_proportions.png"))

# Show the figure
plt.show()
plt.close()

print("Proportions plot saved.")

# Save after computing proportions
adata_file = os.path.join(results_dir, "GFP_combined_proportions_adata.h5ad")
adata.write(adata_file)
print(f"AnnData object saved after computing proportions to {adata_file}")


#############################################################################################################################
# Step 2: Compute RNA velocity
#############################################################################################################################
# Reload the saved AnnData object to continue the analysis
adata = sc.read(os.path.join(results_dir, "GFP_combined_proportions_adata.h5ad"))

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

# Save after computing moments
adata_file = os.path.join(results_dir, "GFP_combined_moments_adata.h5ad")
adata.write(adata_file)
print(f"AnnData object saved after computing moments to {adata_file}")


# Reload the saved AnnData object to continue the analysis
adata = sc.read(os.path.join(results_dir, "mCherry_combined_moments_adata.h5ad"))
print(f"AnnData object reloaded from {adata}")

# Compute RNA velocity
print("Computing RNA velocity.")
scv.tl.velocity(adata)

# Compute RNA velocity graph
scv.tl.velocity_graph(adata)

# Run scVelo Analysis
print("Running scVelo analysis (recover dynamics).")
scv.tl.recover_dynamics(adata, n_jobs=8) # Recover dynamics

# Save after computing RNA velocity in dynamical mode
dynamical_adata_file = os.path.join(results_dir, "mCherry_combined_dynamical_adata_with_velocity.h5ad")
adata.write(dynamical_adata_file)
print(f"AnnData object saved after computing RNA velocity (dynamical mode) to {dynamical_adata_file}")


# Reload the saved AnnData object to continue the analysis
adata = sc.read(os.path.join(results_dir, "mCherry_combined_dynamical_adata_with_velocity.h5ad"))
print(f"AnnData object reloaded from {adata}")

# Compute RNA velocity (dynamical mode)
print("Computing RNA velocity (dynamical mode).") 
scv.tl.velocity(adata, mode="dynamical") 

# Save after computing RNA velocity
intermediate_adata_file = os.path.join(results_dir, "mCherry_combined_intermediate_adata_with_velocity.h5ad")
adata.write(intermediate_adata_file)
print(f"AnnData object saved after computing RNA velocity to {intermediate_adata_file}")


# Reload the saved AnnData object to continue the analysis
adata = sc.read(os.path.join(results_dir, "GFP_combined_intermediate_adata_with_velocity.h5ad"))
print(f"AnnData object reloaded from {adata}")

# Recompute velocity graph
print("Recomputing velocity graph.") 
scv.tl.velocity_graph(adata)

# Plot RNA velocity
print("Plotting RNA velocity.")

# Check and list embeddings in the AnnData object
print("Available embeddings in the AnnData object:")
print(adata.obsm.keys())

# Check available layers
print(adata.layers)

# Plot RNA velocity with streamlines and a legend on the right
fig, ax = plt.subplots(figsize=(10, 8))  # Adjust the figure size as needed
scv.pl.velocity_embedding_stream(adata, basis='umap', color='CellType', ax=ax, legend_loc='right', legend_fontsize=12)
#legend_loc='right', legend_fontsize=12

# Adjust the legend position
legend = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

# Add a title
ax.set_title('RNA Velocity Streamline Plot for GFP', fontsize=16)

# Adjust layout
plt.tight_layout()

# Save the figure
plt.savefig(os.path.join(results_dir, 'GFP_combined_velocity_streamline_plot.png'), bbox_inches='tight')

# Show the figure
plt.show()
plt.close()

print("Streamline velocity plot saved.")

# Save the AnnData object
adata_file = os.path.join(results_dir, "GFP_combined_adata_with_velocity.h5ad")
adata.write(adata_file)
print(f"AnnData object saved to {adata_file}")

