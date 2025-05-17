# Get Special antigen sequence from skin 
python Antigen_Sequence_Selection/src/Get_Special_antigen_data.py --help

# Acquire antigen sequence according to structure of skin
nohup \ 
python Antigen_Sequence_Selection/src/Get_Special_antigen_data.py fetch-fasta \
    --query "(reviewed:true) AND (organism_id:9606) AND (cc_tissue_specificity:("skin" OR "epidermis" OR "dermis" OR "hypodermis")) " \
    --output-fasta "Antigen_Sequence_Selection/antigen_sequnce/Skin/antigen_proteins.fasta" \
    > Antigen_Sequence_Selection/src/run.log 2>&1 &

# Acquire the metadata of antigen sequence according to structure of skin
nohup \
python Antigen_Sequence_Selection/src/Get_Special_antigen_data.py filter-data \
    --input-tsv_url "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Ccc_tissue_specificity%2Ccc_subcellular_location%2Cft_transmem&format=tsv&query=%28%28reviewed%3Atrue%29+AND+%28organism_id%3A9606%29+AND+%28cc_tissue_specificity%3A%28%22skin%22+OR+%22epidermis%22+OR+%22dermis%22+OR+%22hypodermis%22%29%29%29"\
    --filter-pattern "Secreted|Membrane" \
    --output-tsv "Antigen_Sequence_Selection/antigen_sequnce/Skin/proteins_meta.tsv" \
    --filted_output_tsv "Antigen_Sequence_Selection/antigen_sequnce/Skin/filtered_proteins_meta.tsv" \
    >> Antigen_Sequence_Selection/src/run.log 2>&1 &

# Acquire antigen sequence according to cell types
nohup \ 
python Antigen_Sequence_Selection/src/Get_Special_antigen_data.py fetch-fasta \
    --query '(reviewed:true) AND (organism_id:9606) AND (cc_tissue_specificity:("keratinocyte" OR "melanocyte" OR "langerhans cell" OR "merkel cell" OR "fibroblast" OR "hair follicle" OR "sebaceous gland" OR  "sweat gland" OR "mast cell" OR "macrophage" OR "dendritic cell" OR "endothelial cell" OR "pericyte")) ' \
    --output-fasta "Antigen_Sequence_Selection/antigen_sequnce/Blood/antigen_proteins.fasta" \
    >> Antigen_Sequence_Selection/src/run.log 2>&1 &

# Acquire the metadata of antigen sequence according to cell types
nohup \
python Antigen_Sequence_Selection/src/Get_Special_antigen_data.py filter-data \
    --input-tsv_url "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Ccc_tissue_specificity%2Cft_transmem%2Ccc_subcellular_location&format=tsv&query=%28%28reviewed%3Atrue%29+AND+%28organism_id%3A9606%29+AND+%28cc_tissue_specificity%3A%28%22keratinocyte%22+OR+%22melanocyte%22+OR+%22langerhans+cell%22+OR+%22merkel+cell%22+OR+%22fibroblast%22+OR+%22hair+follicle%22+OR+%22sebaceous+gland%22+OR+%22sweat+gland%22+OR+%22mast+cell%22+OR+%22macrophage%22+OR+%22dendritic+cell%22+OR+%22endothelial+cell%22+OR+%22pericyte%22%29%29%29"\
    --filter-pattern "Secreted|Membrane" \
    --output-tsv "Antigen_Sequence_Selection/antigen_sequnce/Blood/proteins_meta.tsv" \
    --filted_output_tsv "Antigen_Sequence_Selection/antigen_sequnce/Blood/filtered_proteins_meta.tsv" \
    >> Antigen_Sequence_Selection/src/run.log 2>&1 &

# Acquire antigen sequence according to THE HUMAN PROTEIN ATLAS
nohup \ 
python Antigen_Sequence_Selection/src/Get_Special_antigen_data.py fetch-fasta2 \
    --query "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/lRyWnJkcCk?compressed=true&format=fasta&query=%28model_organism%3A9606%29" \
    --output-fasta "Antigen_Sequence_Selection/antigen_sequnce/Skin_atlas/antigen_proteins.fasta" \
    >> Antigen_Sequence_Selection/src/run.log 2>&1 &

# Acquire the metadata of antigen sequence according to THE HUMAN PROTEIN ATLAS
nohup \
python Antigen_Sequence_Selection/src/Get_Special_antigen_data.py filter-data \
    --input-tsv_url "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/lRyWnJkcCk?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Ccc_tissue_specificity%2Cft_transmem%2Ccc_subcellular_location&format=tsv&query=%28model_organism%3A9606%29"\
    --filter-pattern "Secreted|Membrane" \
    --output-tsv "Antigen_Sequence_Selection/antigen_sequnce/Skin_atlas/proteins_meta.tsv" \
    --filted_output_tsv "Antigen_Sequence_Selection/antigen_sequnce/Skin_atlas/filtered_proteins_meta.tsv" \
    >> Antigen_Sequence_Selection/src/run.log 2>&1 &

# Split the fasta file to 50 chunks
# 需要处理的子目录列表
target_dirs=("Blood" "Skin" "Skin_atlas")
base_dir="Antigen_Sequence_Selection/antigen_sequnce"
log_file="Antigen_Sequence_Selection/antigen_sequnce/run.log"
for dir in "${target_dirs[@]}"; do
    input_fasta="${base_dir}/${dir}/antigen_proteins.fasta"
    output_dir="${base_dir}/${dir}/"
    
    nohup python Antigen_Sequence_Selection/src/Get_Special_antigen_data.py split-fasta \
        --input-fasta "$input_fasta" \
        --output-dir "$output_dir" \
        --chunk-size 50 \
        >> "$log_file" 2>&1 &
done

# Extract the predicted epitopes of the antigen sequence
# 需要处理的子目录列表
max_jobs=10  # 并发控制
background_jobs=0
target_dirs=("Blood" "Skin" "Skin_atlas")
sub_target_dirs=(p{1..9})  # 生成 p1 p2 p3 ... p9
base_dir="Antigen_Sequence_Selection/antigen_Epitopes_Prediction/BepiPred3"
base_out_dir="Antigen_Sequence_Selection/antigen_Epitopes_Extraction"
log_file="Antigen_Sequence_Selection/antigen_Epitopes_Extraction/run.log"

for dir in "${target_dirs[@]}"; do
    for sub_dir in "${sub_target_dirs[@]}"; do
        mkdir -p "${base_out_dir}/${dir}/${sub_dir}"  # 自动创建嵌套目录
        input_fasta="${base_dir}/${dir}/${sub_dir}/bepipred3_results/Bcell_epitope_preds.fasta"
        output_dir="${base_out_dir}/${dir}/${sub_dir}/epitopes.tsv"
        (   
            echo "Processing ${base_dir}/${dir}/${sub_dir} at $(date)"
            cmd_args=(
                "process-epitopes"
                "--input-fasta" "${input_fasta}"
                "--min-length" "5"
                "--merge-window" "5"
                "--output-tsv" "${output_dir}"
            )
            python "Antigen_Sequence_Selection/src/Get_Process_epitopes.py" "${cmd_args[@]}" 
        ) >> "$log_file" 2>&1 &
        ((background_jobs++))
        if (( background_jobs % max_jobs == 0 )); then
            wait
        fi            
    done
done
wait  
echo "All tasks completed"
