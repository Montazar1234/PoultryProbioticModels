#Code editted from "https://colab.research.google.com/github/Gibbons-Lab/isb_course_2023/blob/main/micom_2023.ipynb"
import pandas as pd
import re
import micom
import os
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import set_start_method
from micom.workflows import build, grow, save_results
from micom.qiime_formats import load_qiime_medium
from micom.workflows import tradeoff
from micom.viz import plot_exchanges_per_taxon
from micom.viz import plot_growth, plot_exchanges_per_sample
from micom.viz import plot_tradeoff


# Set the start method for multiprocessing
set_start_method('spawn', force=True)


# Define the probiotic organism once here
PROBIOTIC_ORGANISM = 'Penicillium_Pd1_GCF_000315645.1'
PROBIOTIC_REL_ABUNDANCE = 0.05


def reshape_microbiome_data(file_path):
    # Load the data
    df = pd.read_excel(file_path, engine='odf')
   
    # Melt the dataframe to long format
    df_melted = df.melt(id_vars=['genus'], var_name='sample_metric', value_name='abundance')
   
    # Split the sample_metric into sample_id and metric (in this case, it will be 'Mean' for all)
    df_melted[['sample_id', 'metric']] = df_melted['sample_metric'].str.rsplit('_', n=1, expand=True)
   
    # Drop NaN values and rows where value is 0
    df_melted = df_melted.dropna().loc[df_melted['abundance'] != 0].reset_index(drop=True)
   
    # Extract data where metric is 'Mean' (this step may be redundant as all metrics are 'Mean')
    df_mean = df_melted.loc[df_melted['metric'] == 'Mean'].copy()
   
    # Add an id column
    df_mean['id'] = df_mean.index + 1
   
    # Rename and reorder columns
    df_final = df_mean.loc[:, ['id', 'sample_id', 'abundance', 'genus']]
    df_final['relative'] = df_final['abundance']
   
    return df_final

def introduce_probiotic(data, probiotic, probiotic_rel):
    """
    Introduce a probiotic organism into the microbial community.
    """
    introduced = pd.DataFrame()
    for smp, df in data.groupby(by='sample_id'):
        df = df[df.relative > 0.00001].copy()
        df = df[df.genus != PROBIOTIC_ORGANISM]
        abund = df.abundance.sum() * probiotic_rel / (1 - probiotic_rel)

        info = df.iloc[0, :].copy()
        info.genus = probiotic
        info.id = probiotic
        info.abundance = abund
        df.loc[df.shape[0] + 1, info.index] = info.values
        df.relative = df.abundance.apply(lambda x: x / df.abundance.sum())
        introduced = pd.concat([introduced, df])
    return introduced

def remove_probiotic(data, probiotic):
    """
    Remove the probiotic organism from the microbial community.
    """
    removed_probiotic = pd.DataFrame()
    for smp, df in data.groupby(by='sample_id'):
        df = df[df.relative > 0.00001].copy()
        df = df[df.genus != probiotic]

        # Re-calculate relative abundances
        df.relative = df.abundance.apply(lambda x: x / df.abundance.sum())

        # Keep consistent with `introduce_probiotic` function
        df['id'] = df['genus']

        removed_probiotic = pd.concat([removed_probiotic, df])
    return removed_probiotic

def prepare_for_building(data):
    # Drop the 'id' column if it exists
    if 'id' in data.columns:
        data = data.drop('id', axis=1)
   
    # Ensure that all object-type columns are strings and replace NaNs
    for col in data.select_dtypes(include=['object']).columns:
        data[col] = data[col].fillna('Unknown').astype(str)
   
    return data


if __name__ == "__main__":
    # Load the manifest and metagenomic data
    manifest = pd.read_csv("./models/manifest.csv")
    metagenomic_file_path = "./aggregated_relative_abundances.ods"

    # Filter for Fungi kingdom
    fungi_df = manifest[manifest['kingdom'] == "Fungi"]

    # Reshape the metagenomic data
    reshaped_data = reshape_microbiome_data(metagenomic_file_path)
   
    medium = pd.read_csv("completed_cornsoymix_new.csv")
    medium = medium.set_index("metabolite")

    # Loop through each fungal organism in the manifest
    for _, fungus in fungi_df.iterrows():
        try:
            current_probiotic = fungus['genus']

            # Skip if the probiotic is already the one defined in the script
            if current_probiotic == PROBIOTIC_ORGANISM:
                continue

            # Define the path for the output zip file
            zip_file_path = f"./cornsoymix_{current_probiotic}/growth_with_{current_probiotic}.zip"

            # Check if the zip file already exists, if so, skip to the next probiotic
            if os.path.exists(zip_file_path):
                print(f"Zip file for {current_probiotic} already exists. Skipping...")
                continue

            # If the zip file doesn't exist, proceed with probiotic introduction and community building
            introduced_data = introduce_probiotic(reshaped_data, current_probiotic, PROBIOTIC_REL_ABUNDANCE)
            introduced_data_prepared = prepare_for_building(introduced_data)

            # Build and grow the microbial community with the current fungus as probiotic
            out_folder_introduced = f"cornsoymix_{current_probiotic}"
            manifest_introduced = build(introduced_data_prepared, model_db="models", out_folder=out_folder_introduced, solver="cplex", cutoff=0.00001, threads=32)
            growth_results = grow(manifest_introduced, out_folder_introduced, medium, tradeoff=0.7, threads=32)

            # Save the growth results
            save_results(growth_results, zip_file_path)
           
        except Exception as e:
            print(f"An error occurred with {current_probiotic}: {e}")
            continue

