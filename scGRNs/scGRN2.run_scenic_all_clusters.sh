#!/bin/bash

# ----------------------------
# Activate conda environment
# ----------------------------
source ~/.bashrc
conda activate pyscenic

# ----------------------------
# Constants
# ----------------------------
# Using the requested base directory.
BASE_DIR="/home/greenhamlab/snTC/20240411.GRNs.ParallelShuffled"
DB_DIR="/home/greenhamlab/scPlant/scPlantData/cisTarget_databases"
TF_FILE="${DB_DIR}/At_tfs.CyclingAdded.txt"
RANKING_FILE="${DB_DIR}/At.genes_vs_motifs.rankings.feather"
MOTIF_ANNOT_FILE="${DB_DIR}/At_motif2TF.tbl"
NUM_WORKERS=32
SEED=777

# ----------------------------
# Cluster list (excluding 16)
# ----------------------------
clusters=("0" "1" "5" "2" "11" "3" "4" "12_2" "10" "9" "14" "6" "12_0" "7" "8_0" "12_1" "13" "8_1" "15")

# ----------------------------
# Loop through clusters and process loom files
# ----------------------------
for cluster in "${clusters[@]}"; do
    echo "Processing cluster: $cluster"
    CLUSTER_DIR="${BASE_DIR}/GRNs.c${cluster}.ClusterOnlyCycling.c${cluster}OnlyCells"
    
    if [ ! -d "$CLUSTER_DIR" ]; then
        echo "Error: Cluster directory ${CLUSTER_DIR} does not exist." >&2
        continue
    fi
    
    # For debugging: list files in the cluster directory
    echo "Listing files in ${CLUSTER_DIR}:"
    ls -l "$CLUSTER_DIR"

    # Build an array of all loom files (original plus shuffled) in the cluster directory
    loom_files=("$CLUSTER_DIR"/exprMat*.loom)
    
    # Check that at least one file was found
    if [ ! -e "${loom_files[0]}" ]; then
        echo "Error: No loom files found in ${CLUSTER_DIR} matching exprMat*.loom" >&2
        continue
    fi

    # Iterate over each loom file in the array
    for loom_file in "${loom_files[@]}"; do
        loom_basename=$(basename "$loom_file" .loom)
        echo "  Processing loom file: $loom_basename"

        # Create an output subdirectory for the current loom file's scenic output
        OUTPUT_DIR="${CLUSTER_DIR}/${loom_basename}_scenic"
        mkdir -p "$OUTPUT_DIR"

        # Define output file names
        ADJ_OUT="${OUTPUT_DIR}/adj.tsv"
        REG_OUT="${OUTPUT_DIR}/reg.tsv"
        LOOM_OUT="${OUTPUT_DIR}/pyscenicOutput.loom"

        # Step 1: arboreto
        echo "    Running arboreto for loom file $loom_basename"
        arboreto_with_multiprocessing.py \
            "$loom_file" \
            "$TF_FILE" \
            --method genie3 \
            --output "$ADJ_OUT" \
            --num_workers $NUM_WORKERS \
            --seed $SEED

        # Step 2: ctx
        echo "    Running ctx for loom file $loom_basename"
        pyscenic ctx \
            "$ADJ_OUT" \
            "$RANKING_FILE" \
            --annotations_fname "$MOTIF_ANNOT_FILE" \
            --expression_mtx_fname "$loom_file" \
            --no_pruning \
            --output "$REG_OUT" \
            --num_workers $NUM_WORKERS

        # Step 3: aucell
        echo "    Running aucell for loom file $loom_basename"
        pyscenic aucell \
            "$loom_file" \
            "$REG_OUT" \
            --output "$LOOM_OUT" \
            --num_workers $NUM_WORKERS \
            --seed $SEED

        echo "  Finished processing loom file: $loom_basename"
        echo "---------------------------------------------"
    done

    echo "Finished cluster $cluster"
    echo "============================================="
done


# ----------------------------
# Constants
# ----------------------------
# Using the requested base directory.
BASE_DIR="/home/greenhamlab/snTC/20240411.GRNs.ParallelShuffled"
DB_DIR="/home/greenhamlab/scPlant/scPlantData/cisTarget_databases"
TF_FILE="${DB_DIR}/At_tfs.CyclingAdded.txt"
RANKING_FILE="${DB_DIR}/At.genes_vs_motifs.rankings.feather"
MOTIF_ANNOT_FILE="${DB_DIR}/At_motif2TF.tbl"
NUM_WORKERS=32
SEED=777

# ----------------------------
# Cluster list (only missed sub‑clusters)
# ----------------------------
clusters=("8_0" "8_1" "12_0" "12_1" "12_2")

# ----------------------------
# Loop through clusters and process loom files
# ----------------------------
for cluster in "${clusters[@]}"; do
    echo "Processing cluster: $cluster"
    CLUSTER_DIR="${BASE_DIR}/GRNs.c${cluster}.ClusterOnlyCycling.c${cluster}OnlyCells"
    
    if [ ! -d "$CLUSTER_DIR" ]; then
        echo "Error: Cluster directory ${CLUSTER_DIR} does not exist." >&2
        continue
    fi
    
    # For debugging: list files in the cluster directory
    echo "Listing files in ${CLUSTER_DIR}:"
    ls -l "$CLUSTER_DIR"

    # Build an array of all loom files (original plus shuffled) in the cluster directory
    loom_files=("$CLUSTER_DIR"/exprMat*.loom)
    
    # Check that at least one file was found
    if [ ! -e "${loom_files[0]}" ]; then
        echo "Error: No loom files found in ${CLUSTER_DIR} matching exprMat*.loom" >&2
        continue
    fi

    # Iterate over each loom file in the array
    for loom_file in "${loom_files[@]}"; do
        loom_basename=$(basename "$loom_file" .loom)
        echo "  Processing loom file: $loom_basename"

        # Create an output subdirectory for the current loom file's scenic output
        OUTPUT_DIR="${CLUSTER_DIR}/${loom_basename}_scenic"
        mkdir -p "$OUTPUT_DIR"

        # Define output file names
        ADJ_OUT="${OUTPUT_DIR}/adj.tsv"
        REG_OUT="${OUTPUT_DIR}/reg.tsv"
        LOOM_OUT="${OUTPUT_DIR}/pyscenicOutput.loom"

        # Step 1: arboreto
        echo "    Running arboreto for loom file $loom_basename"
        arboreto_with_multiprocessing.py \
            "$loom_file" \
            "$TF_FILE" \
            --method genie3 \
            --output "$ADJ_OUT" \
            --num_workers $NUM_WORKERS \
            --seed $SEED

        # Step 2: ctx
        echo "    Running ctx for loom file $loom_basename"
        pyscenic ctx \
            "$ADJ_OUT" \
            "$RANKING_FILE" \
            --annotations_fname "$MOTIF_ANNOT_FILE" \
            --expression_mtx_fname "$loom_file" \
            --no_pruning \
            --output "$REG_OUT" \
            --num_workers $NUM_WORKERS

        # Step 3: aucell
        echo "    Running aucell for loom file $loom_basename"
        pyscenic aucell \
            "$loom_file" \
            "$REG_OUT" \
            --output "$LOOM_OUT" \
            --num_workers $NUM_WORKERS \
            --seed $SEED

        echo "  Finished processing loom file: $loom_basename"
        echo "---------------------------------------------"
    done

    echo "Finished cluster $cluster"
    echo "============================================="
done

echo "All clusters processed."

