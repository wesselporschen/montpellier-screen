#!/bin/bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate ringtail

for db in *.db; do 
    name="${db%.db}"
    echo "Processing $name"

    rt_process_vs read --input_db "${name}.db" --ligand_name Rahimova Ghoteimi Grossjean --bookmark_name "all_actives_${name}"
    rt_process_vs read --input_db "${name}.db" --ligand_name ZINC --bookmark_name "all_decoys_${name}"

    rt_process_vs read --input_db "${name}.db" -xs "all_actives_${name}"
    rt_process_vs read --input_db "${name}.db" -xs "all_decoys_${name}"

done 



