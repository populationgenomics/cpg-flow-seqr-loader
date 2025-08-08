#!/usr/bin/env bash

### I heard you like scripts, so I put a script in your script so you can run a script while you run a script

### This should absolutely not be run in main

# take 3 String values as arguments, there must be a clean way to do this
if [ -z "$1" ]; then
  echo "Usage: $0 <input_file> <plot_output_folder> <api_key>"
  exit 1
fi
if [ -z "$2" ]; then
  echo "Usage: $0 <input_file> <plot_output_folder> <api_key>"
  exit 1
fi
if [ -z "$3" ]; then
  echo "Usage: $0 <input_file> <plot_output_folder> <api_key>"
  exit 1
fi

INPUT_FILE="$1"
OUTPUT_DIR="$2"
API_KEY="$3"

# stage outputs locally, then copy out
mkdir alphagenome_outputs

python3 src/cpg_seqr_loader/scripts/mw_edit.py \
    --variants "${INPUT_FILE}" \
    --output alphagenome_outputs \
    --api-key "${API_KEY}"

gcloud storage cp -r alphagenome_outputs "${OUTPUT_DIR}"
