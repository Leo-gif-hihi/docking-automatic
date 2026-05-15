#!/bin/bash

# Arguments
PROTEIN_FILE=$1
LIGAND_FILE=$2
BOX_FILE=$3
OUTPUT_DIR=$4
CPUS=${5:-0}

# Extract paths and basenames
PROTEIN_DIR=$(dirname "$PROTEIN_FILE")
LIGAND_DIR=$(dirname "$LIGAND_FILE")

PROTEIN_BASE=$(basename "$PROTEIN_FILE")
PROTEIN_BASE="${PROTEIN_BASE%.*}"

LIGAND_BASE=$(basename "$LIGAND_FILE")
LIGAND_BASE="${LIGAND_BASE%.*}"

# Intermediate and Output paths
# mk_prepare_receptor appends .pdbqt to the -o value, so we strip it.
PROTEIN_PREP_OUT="${PROTEIN_DIR}/${PROTEIN_BASE}" 
PROTEIN_PDBQT="${PROTEIN_DIR}/${PROTEIN_BASE}.pdbqt"
LIGAND_PDBQT="${LIGAND_DIR}/${LIGAND_BASE}.pdbqt"

OUT_PDBQT="${OUTPUT_DIR}/${PROTEIN_BASE}_${LIGAND_BASE}_vina_out.pdbqt"
OUT_SDF="${OUTPUT_DIR}/${PROTEIN_BASE}_${LIGAND_BASE}_vina_out.sdf"
OUT_LOG="${OUTPUT_DIR}/${PROTEIN_BASE}_${LIGAND_BASE}_vina.log"

echo "Preparing Receptor..."
mk_prepare_receptor.py -i "$PROTEIN_FILE" -o "$PROTEIN_PREP_OUT" -p

echo "Preparing Ligand..."
mk_prepare_ligand.py -i "$LIGAND_FILE" -o "$LIGAND_PDBQT"

VINA_CPU_ARG=""
if [ "$CPUS" -gt 0 ]; then
    VINA_CPU_ARG="--cpu $CPUS"
fi

echo "Running Vina..."
vina --receptor "$PROTEIN_PDBQT" --ligand "$LIGAND_PDBQT" \
       --config "$BOX_FILE" \
       --exhaustiveness=32 $VINA_CPU_ARG --out "$OUT_PDBQT" > "$OUT_LOG"

echo "Exporting to SDF..."
mk_export.py "$OUT_PDBQT" -s "$OUT_SDF"