import argparse
from rdkit import Chem
import pandas as pd 
import os



def extract_actives_inactives_sdf():
    parser = argparse.ArgumentParser()
    parser.add_argument("labels_csv")
    parser.add_argument("sdf_file")

    args = parser.parse_args()

    labels = pd.read_csv(args.labels_csv)
    labels_dict = dict(zip(labels["ID"], labels["ACTIVE_LABEL"]))

    # Input / output files

    supplier = Chem.SDMolSupplier(args.sdf_file)
    active_writer = Chem.SDWriter("true_actives.sdf")
    inactive_writer = Chem.SDWriter("true_inactives.sdf")

    n_active = 0
    n_inactive = 0
    n_missing = 0

    for mol in supplier:
        if mol is None:
            continue

        lig_id = mol.GetProp("_Name")

        if lig_id not in labels_dict:
            n_missing += 1
            continue

        if labels_dict[lig_id]:
            active_writer.write(mol)
            n_active += 1

        else:
            inactive_writer.write(mol)
            n_inactive += 1

    active_writer.close()
    inactive_writer.close()

    print(f"Actives: {n_active}")
    print(f"Inactives: {n_inactive}")
    print(f"Missing labels: {n_missing}")
