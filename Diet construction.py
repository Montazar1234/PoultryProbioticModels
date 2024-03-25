import pandas as pd
from micom.workflows.db_media import complete_db_medium, check_db_medium

def main():
    # Read the diet file and set the index to "metabolite"
    diet = pd.read_csv("./Example start diet.csv")
    diet.set_index("metabolite", inplace=True)
    
    # Read the annotations file
    annotations = pd.read_csv("./annotations.csv")
    
    # Reset the index and create a reaction column
    diet.reset_index(inplace=True)
    diet["reaction"] = "EX_" + diet.metabolite + "(e)"
    
    # Merge the diet and annotations dataframes
    skeleton = pd.merge(diet, annotations, on="metabolite")
    
    # Rename columns
    skeleton.rename(columns={"reaction_x": "reaction"}, inplace=True)
    
    # Update the global_id and reaction columns
    skeleton["global_id"] = skeleton.reaction
    skeleton["reaction"] = "EX_" + skeleton.metabolite + "_m"
    
    # Complete the medium using MiCOM
    manifest, imports = complete_db_medium(
        "models_poultry", 
        medium=skeleton, 
        growth=0.01, 
        threads=3, 
        max_added_import=10,
        strict=["EX_o2(e)"],  
        weights="mass"        
    )
    
    # Calculate added flux
    filled = imports.max()
    added = filled.sum() - skeleton.loc[skeleton.reaction.isin(filled.index), "flux"].sum()
    print(f"Added flux is {added:.2f}/{filled.sum():.2f} mmol/h.")
    
    # Identify added fluxes
    added_fluxes = filled.copy()
    shared = added_fluxes.index[added_fluxes.index.isin(skeleton.reaction)]
    added_fluxes[shared] -= skeleton.flux[shared]
    
    # Create completed medium dataframe
    added_df = filled[filled > 1e-8].reset_index()
    added_df.iloc[:, 0] = added_df.iloc[:, 0].str.replace("EX_|_m$", "", regex=True)
    added_df.columns = ["metabolite", "flux"]
    completed = pd.merge(added_df, annotations, on="metabolite", how="left")
    completed["reaction"] = "EX_" + completed.metabolite + "_m"
    completed["global_id"] = "EX_" + completed.metabolite + "(e)"
    
    # Check the completed medium
    check = check_db_medium("models_poultry", medium=completed, threads=3)
    print(check.growth_rate.describe())
    
    # Save the completed medium to a CSV file
    completed.to_csv("Example gap filled diet.csv", index=False)

if __name__ == "__main__":
    main()
