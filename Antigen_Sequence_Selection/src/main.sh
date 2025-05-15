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

