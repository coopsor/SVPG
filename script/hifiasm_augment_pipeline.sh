#!/bin/bash
# This script runs hifiasm to assemble genomes from FASTA files in directories starting with "HG", so you need to run it in the parent directory containing those directories.

input_dir="./"
output_dir="./"
hifiasm="/home/huheng/download/hifiasm/hifiasm"
gfatools="/home/huheng/download/gfatools/gfatools"

start_time=$(date +%s)

# Loop through all directories starting with "HG"
for dir in HG*; do
    # Enter the directory
    cd "$dir" || continue

    # Use the directory name as the prefix
    prefix=$(basename "$dir")

    # Run hifiasm to assemble
    $hifiasm -o $prefix -t 128 ${prefix}.fasta

    # Convert GFA to FASTA
    $gfatools gfa2fa ${prefix}.bp.p_ctg.gfa > assem.fa

    # Return to the parent directory
    cd ..
done

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Time taken of hifiasm for $dir: $elapsed_time seconds"

# Prepare list of assembly FASTA files
assem_files=()
for dir in HG*; do
    if [[ -f "$dir/assem.fa" ]]; then
        # Prefix sequence headers with directory name
#        sed "s/^>/&${dir}_/" "$dir/assem.fa" > "$dir/assem.fa.tmp" && mv "$dir/assem.fa.tmp" "$dir/assem.fa"
        assem_files+=("$dir/assem.fa")
    fi
done

# Build the pangenome graph using minigraph
minigraph -cxggs -t128 /data/huheng/dataset/human/pangenome/GRCh38-90c.r518.gfa "${assem_files[@]}" > hifiasm.110c_38_new.gfa

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Time of minigraph taken for $dir: $elapsed_time seconds"