import MDAnalysis as mda
import prolif as plf
from rdkit import Chem
import os
from collections import defaultdict
import pandas as pd
import sys
import json
import pickle
from concurrent.futures import ProcessPoolExecutor

# Collect all sdf files for actives, inactives, and decoys in nested dictionary
results_dir = "./results"
sdf_data = {}

for root, dirs, files in os.walk(results_dir):
    rel_path = os.path.relpath(root, results_dir)

    current = sdf_data

    if "_sdf" in rel_path and '10poses' not in rel_path:
        for part in rel_path.split(os.sep):
            current = current.setdefault(part, {})

        for file in files:
            if file.endswith(".sdf"):
                current[file] = [os.path.join(root, file), None] # Store the path and a placeholder for the DataFrame

protein_conformers = {
    '4h2i_c_1_dimer_sdf': './results/12Jun_rg/structures_Rg_box/4h2i_c_1/4h2i_c_1_mae.pdb',
    '4h2i_c_2_dimer_sdf': './results/12Jun_rg/structures_Rg_box/4h2i_c_2/4h2i_c_2_mae.pdb',
    '4h2i_c_3_dimer_sdf': './results/12Jun_rg/structures_Rg_box/4h2i_c_3/4h2i_c_3_mae.pdb',
    '4h2i_c_4_dimer_sdf': './results/12Jun_rg/structures_Rg_box/4h2i_c_4/4h2i_c_4_mae.pdb',
    '4h2i_crys_dimer_sdf': './results/12Jun_rg/structures_Rg_box/4h2i_crys/4h2i_crys_mae.pdb',
    '6tve_c_1_dimer_sdf': './results/12Jun_rg/structures_Rg_box/6tve_c_1/6tve_c_1_mae.pdb',
    '6tve_c_2_dimer_sdf': './results/12Jun_rg/structures_Rg_box/6tve_c_2/6tve_c_2_mae.pdb',
    '6tve_c_3_dimer_sdf': './results/12Jun_rg/structures_Rg_box/6tve_c_3/6tve_c_3_mae.pdb',
    '6tve_c_4_dimer_sdf': './results/12Jun_rg/structures_Rg_box/6tve_c_4/6tve_c_4_mae.pdb',
    '6tve_crys_dimer_sdf': './results/12Jun_rg/structures_Rg_box/6tve_crys/6tve_crys_mae.pdb'
}

# Generate pandas dataframe of interaction fingerprints for sdf file and protein conformer
def process_sdf_files(ligand_sdf, protein_conformer) -> pd.DataFrame:
    fp_ligands = plf.sdf_supplier(ligand_sdf)
    protein = plf.Molecule(Chem.MolFromPDBFile(protein_conformer))

    interaction_fingerprint = plf.Fingerprint()
    interaction_fingerprint.run_from_iterable(fp_ligands, protein)

    df = interaction_fingerprint.to_dataframe()
    return df




processed = {
    "12Jun_g543",
    "12Jun_rg",
    "12Jun_rg10poses",
    "16Jun_rg2",
    "16Jun_rg3",
    "5Jun_dimer",
}


def worker(task):
    conformer, ligand_sdf_key, ligand_sdf, protein_conformer = task

    print(f"Generating fingerprint dataframe for {conformer}: {ligand_sdf_key}")

    df = process_sdf_files(
        ligand_sdf=ligand_sdf,
        protein_conformer=protein_conformer,
    )

    return conformer, ligand_sdf_key, df


if __name__ == "__main__":

    for docking_run_key, conformers in sdf_data.items():

        if docking_run_key in processed:
            print(f"Already processed docking run {docking_run_key}. Skipping.")
            continue

        print(f"\nProcessing {docking_run_key}")

        # Build the task list for this docking run
        tasks = []

        for conformer, ligand_sdfs in conformers.items():
            protein_conformer = protein_conformers[conformer]

            for ligand_sdf_key, data in ligand_sdfs.items():
                tasks.append(
                    (
                        conformer,
                        ligand_sdf_key,
                        data[0],              # ligand sdf filename
                        protein_conformer,
                    )
                )

        # Process all SDFs in parallel
        with ProcessPoolExecutor() as executor:
            for conformer, ligand_sdf_key, df in executor.map(worker, tasks):
                conformers[conformer][ligand_sdf_key][1] = df

        # Save this docking run immediately
        with open(f"./prolif_fps/{docking_run_key}.pkl", "wb") as f:
            pickle.dump({docking_run_key: conformers}, f)

        print(f"Saved {docking_run_key}.pkl")
