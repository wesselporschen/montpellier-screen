import argparse
from rdkit import Chem
import pandas as pd 
import os

def extract_actives_inactives_sdf(labels_csv: str, 
                                  sdf_file: str,
                                  active_sdf: str,
                                  inactive_sdf: str) -> None:

    print(f"Extracting true active and inactive sdfs")
    print(f"{labels_csv=}")
    print(f"{sdf_file=}")
    print(f"{active_sdf=}")
    print(f"{inactive_sdf=}")

    labels = pd.read_csv(labels_csv)
    labels_dict = dict(zip(labels["ID"], labels["ACTIVE_LABEL"]))

    # Input / output files

    supplier = Chem.SDMolSupplier(sdf_file)
    active_writer = Chem.SDWriter(active_sdf)
    inactive_writer = Chem.SDWriter(inactive_sdf)

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



# CONFIG ------

CWD = os.getcwd()
dirs = os.listdir()
# Change to dockdirs that need to be done, otherwise every docking run dir gets processed again
#dockdirs = [dir for dir in dirs if "Jun" in dir]
dockdirs = ['23Jun_g543_3', '23Jun_g543_2']
print(f"{dockdirs=}")

for dockdir in dockdirs:
    os.chdir(f"{CWD}/{dockdir}")

    subdirs = os.listdir()
    conformerdirs = [dir for dir in subdirs if "_dimer_sdf" in dir]

    for conformerdir in conformerdirs:
        os.chdir(f"{CWD}/{dockdir}/{conformerdir}")

        print(f"Directory: {dockdir} {conformerdir}")

        # find montpellier sdf 
        montpellier_matches = [f for f in os.listdir(".") if f.endswith("montpellier.sdf")]
        montpellier_sdf = montpellier_matches[0]

        # set active and inactive sdf outputs 
        conformerbase = conformerdir.replace("_sdf", "")
        active_sdf = f"{conformerbase}_true_actives.sdf"
        inactive_sdf = f"{conformerbase}_true_inactives.sdf"

        #print(montpellier_sdf, active_sdf, inactive_sdf)

        extract_actives_inactives_sdf(labels_csv='/Users/wes/Documents/bps/RP2/cd73_dualantagonist/docking/montpellier-screen/results/montpellier_ligands.csv',
                                      sdf_file=montpellier_sdf,
                                      active_sdf=active_sdf,
                                      inactive_sdf=inactive_sdf)



