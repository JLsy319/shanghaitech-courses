#!/bin/bash

GPU_ID=$1
if [ -z "$GPU_ID" ]; then
  echo "Error：Input GPU ID (0, 1, 3) as the first parameters"
  echo "Usage: ./run_q3_predictions.sh <GPU_ID>"
  exit 1
fi

LOCK_DIR="/tmp/gpu${GPU_ID}.lock"
MODEL_PARAMS_DIR="/home/guest/compbio-2025/simpified_af3/para"
DATABASES_DIR="/home/guest/compbio-2025/simpified_af3/fake_db"
inp="/home/guest/compbio-2025/proj4_q3_precompute"
outp="/home/guest/compbio-2025/lsy/output"
mkdir -p "$outp"

AF3_CONTAINER_SIF="/home/guest/compbio-2025/simpified_af3/alphafold3_20241205.sif"
CONTAINER_INPUT_DIR="/root/af_input"
CONTAINER_OUTPUT_DIR="/root/af_output"
CONTAINER_MODELS_DIR="/root/models"
CONTAINER_DATABASES_DIR="/root/public_databases"

if mkdir "$LOCK_DIR" 2>/dev/null; then
    echo "Successfully lock GPU ${GPU_ID}"
    trap 'rm -rf "$LOCK_DIR"; echo "GPU ${GPU_ID} 已释放"; exit' INT TERM EXIT HUP QUIT
    shopt -s nullglob
    echo "Start dealing with precomputed MSA json files on GPU ${GPU_ID} ..."
    echo "Input directory (precomputed MSA json files): $inp"
    echo "Output directory: $outp"
    echo "Use model parameters: $MODEL_PARAMS_DIR"
    echo "Use database (fake_db): $DATABASES_DIR"
    target_json_files=(
        "or22c-ms/or22c-ms_data.json"
        "or24a-2a/or24a-2a_data.json"
        "or24a-mps/or24a-mps_data.json"
        "or85c-2h/or85c-2h_data.json"
        "or85c-3r/or85c-3r_data.json"
        "or85c-3s/or85c-3s_data.json"
    )
    for filename_to_process in "${target_json_files[@]}"; do
        json_file_path="$inp/$filename_to_process"
        if [ -f "$json_file_path" ]; then
            echo "----------------------------------------"
            echo "Procssing: $filename_to_process"
            echo "----------------------------------------"
            CUDA_VISIBLE_DEVICES=${GPU_ID} apptainer exec \
                --nv \
                -B ${inp}:${CONTAINER_INPUT_DIR} \
                -B ${outp}:${CONTAINER_OUTPUT_DIR} \
                -B ${MODEL_PARAMS_DIR}:${CONTAINER_MODELS_DIR} \
                -B ${DATABASES_DIR}:${CONTAINER_DATABASES_DIR} \
                ${AF3_CONTAINER_SIF} \
                python /app/alphafold/run_alphafold.py \
                    --json_path=${CONTAINER_INPUT_DIR}/${filename_to_process} \
                    --model_dir=${CONTAINER_MODELS_DIR} \
                    --db_dir=${CONTAINER_DATABASES_DIR} \
                    --output_dir=${CONTAINER_OUTPUT_DIR} \
                    --max_template_date="9999-01-01" \
                    --jackhmmer_n_cpu=8 # [cite: 568] (Not important when MSA has been precomputed)
            if [ $? -eq 0 ]; then
                echo "File $filename_to_process is successfully processed."
            else
                echo "Warning: File $filename_to_process is wrongly processed."
            fi
            echo "----------------------------------------"
        else
            echo "Warning: File $json_file_path isn't found. Skip."
        fi
    done

    shopt -u nullglob
    echo "All files are processed."

else
    echo "Error: GPU ${GPU_ID} is occupied."
    exit 1
fi