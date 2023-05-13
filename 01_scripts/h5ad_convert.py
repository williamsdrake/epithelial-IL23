import scanpy as sc
from scipy import io
import sys
import os

h5ad_directory = "/data/williamsdrw/h5ad_convert/"

h5ad_files = [file for file in os.listdir(h5ad_directory) if file.endswith(".h5ad")]

for h5ad_file in h5ad_files:

    # Read the h5ad file
    adata = sc.read_h5ad(os.path.join(h5ad_directory, h5ad_file))

    # Set the output directory
    out_dir = os.path.join(h5ad_directory, h5ad_file.split(".h5ad")[0])

    # Create the output directory
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Write the barcodes to a file called `barcodes.tsv`
    with open(os.path.join(out_dir, "barcodes.tsv"), "w") as f:
        for item in adata.obs_names:
            f.write(item + "\n")

    # Write the gene expression values to a file called `features.tsv`
    with open(os.path.join(out_dir, "features.tsv"), "w") as f:
        for item in ['\t'.join([x,x,'Gene Expression']) for x in adata.var_names]:
            f.write(item + "\n")
    # Generate a matrix file
    io.mmwrite(out_dir +'/matrix', adata.X.T)

    # Write the metadata to a file called `metadata.csv`
    adata.obs.to_csv(os.path.join(out_dir, "metadata.csv"))

    print("File " + h5ad_file + " has been processed")
