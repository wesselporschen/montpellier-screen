#!/bin/bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate ringtail

for db in *.db; do 
    name="${db%.db}"
    echo "Processing $name"

    rt_process_vs read --input_db "${name}.db" --bookmark_name "all_actives_${name}" -sdf "${name}_sdf"
    rt_process_vs read --input_db "${name}.db" --bookmark_name "all_decoys_${name}" -sdf "${name}_sdf"

done 



