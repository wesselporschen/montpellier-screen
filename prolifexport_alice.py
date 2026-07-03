import prolif as plf
from rdkit import Chem
import pandas as pd
from pathlib import Path


conformer_jobs = ["4h2i_c_4_percent5", "6tve_c_1_percent5", "6tve_c_2_percent5"]

protein_structures = {
    '4h2i_c_4_percent5': '/home/s2004267/data1/cd73/docking/vs/structures_Rg/4h2i_c_4/4h2i_c_4_mae.pdb',
    '6tve_c_1_percent5': '/home/s2004267/data1/cd73/docking/vs/structures_Rg/6tve_c_1/6tve_c_1_mae.pdb',
    '6tve_c_2_percent5': '/home/s2004267/data1/cd73/docking/vs/structures_Rg/6tve_c_2/6tve_c_2_mae.pdb',
}

sdf_files = {
    '4h2i_c_4_percent5': '/home/s2004267/data1/cd73/docking/vs/results/dbs/percentile_5/4h2i_c_4_percent5_all_results.sdf',
    '6tve_c_1_percent5': '/home/s2004267/data1/cd73/docking/vs/results/dbs/percentile_5/6tve_c_1_percent5_all_results.sdf',
    '6tve_c_2_percent5': '/home/s2004267/data1/cd73/docking/vs/results/dbs/percentile_5/6tve_c_2_percent5_all_results.sdf',

        }

# Generate pandas dataframe of interaction fingerprints for sdf file and protein conformer
def process_sdf_files(ligand_sdf, protein_conformer) -> pd.DataFrame:
    # Read molecules from the SDF
    rdkit_mols = [mol for mol in Chem.SDMolSupplier(ligand_sdf) if mol is not None]

    # Extract ligand IDs (first line of each SDF record)
    ligand_ids = [
        mol.GetProp("_Name") if mol.HasProp("_Name") else None
        for mol in rdkit_mols
    ]

    fp_ligands = plf.sdf_supplier(ligand_sdf)
    protein = plf.Molecule(Chem.MolFromPDBFile(protein_conformer))

    interaction_fingerprint = plf.Fingerprint()
    interaction_fingerprint.run_from_iterable(fp_ligands, protein)

    df = interaction_fingerprint.to_dataframe()
    df.insert(0, "LigandID", ligand_ids)
    return df


if __name__ == "__main__":

    for conformer_job in conformer_jobs:
        assert Path(sdf_files[conformer_job]).exists()
        assert Path(protein_structures[conformer_job]).exists()

        print(f"Generating prolif df from sdf files for {conformer_job}...")
        prolif_df = process_sdf_files(ligand_sdf=sdf_files[conformer_job],
                                      protein_conformer=protein_structures[conformer_job])

        print("Done.")
        print(f"Saving df as pkl...")
        prolif_df.to_pickle(f"{conformer_job}.pkl")
        print("Done.\n")

