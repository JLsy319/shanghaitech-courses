# Project 4

## Task 1

PurA: P0A7D4  
PurZ: G3FFN6  
PurZ0: A0A7L7SI10

| Protein - Ligand | ipTM | pTM |
| :---: | :---: | :---: |
| PurA-ATP dimer | 0.9 | 0.92 |
| PurA-GTP dimer | 0.91 | 0.94 |
| PurZ-ATP dimer | 0.94 | 0.95 |
| PurZ-GTP dimer | 0.89 | 0.92 |
| PurZ0-ATP dimer | 0.89 | 0.92 |
| PurZ0-GTP dimer | 0.87 | 0.91 |

<img src="pic/q1.png" alt="af3 result" style="zoom:67%;" />

## Task 2

According to the paper, the combinations are:
- ThsB TypeI monomer/dimer + Tad3 monomer/dimer
- ThsB TypeI monomer/dimer + Tad4 monomer/dimer
- ThsB TypeI monomer/dimer + Tad5 monomer/dimerr
- ThsB TypeI monomer/dimerr + Tad6 monomer/dimer
- ThsB TypeII monomer + Tad7 monomer OR ThsB TypeII dimer + Tad7 dimer
- ThsA TypeII monomer + Tad8 dimer
- CD-NTase monomer + Acb3 monomer

### Result of ThsB and Tad3

1. monomer + monomer

<img src="pic/ThsB-and-Tad3-mo-mo.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.5, pTM = 0.69

2. monomer + dimer

<img src="pic/ThsB-and-Tad3-mo-di.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.87, pTM = 0.85

3. dimer + monomer

<img src="pic/ThsB-and-Tad3-di-mo.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.57, pTM = 0.68

4. dimer + dimer

<img src="pic/ThsB-and-Tad3-di-di.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.87, pTM = 0.88

Complex composition: ThsB TypeI dimer + Tad3 dimer

### Result of ThsB and Tad4

1. monomer + monomer

<img src="pic/ThsB-and-Tad4-mo-mo.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.92, pTM = 0.9

2. monomer + dimer

<img src="pic/ThsB-and-Tad4-mo-di.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.57, pTM = 0.66

3. dimer + monomer

<img src="pic/ThsB-and-Tad4-di-mo.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.48, pTM = 0.6

4. dimer + dimer

<img src="pic/ThsB-and-Tad4-di-di.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.44, pTM = 0.54

Complex composition: ThsB TypeI monomer + Tad4 monomer

### Result of ThsB and Tad5

1. monomer + monomer

<img src="pic/ThsB-and-Tad5-mo-mo.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.8, pTM = 0.82

2. monomer + dimer

<img src="pic/ThsB-and-Tad5-mo-di.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.77, pTM = 0.79

3. dimer + monomer

<img src="pic/ThsB-and-Tad5-di-mo.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.51, pTM = 0.63

4. dimer + dimer

<img src="pic/ThsB-and-Tad5-di-di.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.85, pTM = 0.86

Complex composition: ThsB TypeI dimer + Tad5 dimer

### Result of ThsB and Tad6

1. monomer + monomer

<img src="pic/ThsB-and-Tad6-mo-mo.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.83, pTM = 0.79

2. monomer + dimer

<img src="pic/ThsB-and-Tad6-mo-di.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.42, pTM = 0.53

3. dimer + monomer

<img src="pic/ThsB-and-Tad6-di-mo.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.45, pTM = 0.56

4. dimer + dimer

<img src="pic/ThsB-and-Tad6-di-di.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.35, pTM = 0.46

Complex composition: ThsB TypeI monomer + Tad6 monomer

### Result of ThsB and Tad7

1. monomer + monomer

<img src="pic/ThsB-and-Tad7_monomer.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.95, pTＭ = 0.92

2. dimer + dimer

<img src="pic/ThsB-and-Tad7_dimer.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.51, pTM = 0.59

Complex composition: ThsB TypeII monomer + Tad7 monomer

### Result of ThsA and Tad8

<img src="pic/ThsA-and-Tad8.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.81, pTM = 0.82

### Result of CD-NTase and Acb3

<img src="pic/CD-NTase-and-Acb3.jpg" alt="result" style="zoom:50%;" />

ipTM = 0.2, pTM = 0.42

## Task 3

```shell
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
```

### or22c and methyl salicylate

```json
{
    "name": "Or22c_Methyl-salicylate_complex",
    "modelSeeds": [1],
    "sequences": [
      {
        "protein": {
          "id": "A",
          "sequence": "precomputed MSA"
        }
      },
      {
        "ligand": {
          "id": "L",
          "smiles": "COC(=O)C1=CC=CC=C1O"
        }
      }
    ],
    "dialect": "alphafold3",
    "version": 1
  }
```

<img src="pic/or22c-ms.png" alt="result" style="zoom:50%;" />

ipTM = 0.92, pTM = 0.81

### or24a and 2-acetylpyridine

```json
{
    "name": "Or24a_2-acetylpyridine_complex",
    "modelSeeds": [1],
    "sequences": [
      {
        "protein": {
          "id": "A",
          "sequence": "precomputed MSA"
        }
      },
      {
        "ligand": {
          "id": "L",
          "smiles": "CC(=O)C1=CC=CC=N1"
        }
      }
    ],
    "dialect": "alphafold3",
    "version": 1
  }
```

<img src="pic/or24a-2a.png" alt="result" style="zoom:50%;" />

ipTM = 0.76, pTM = 0.7

### or24a and methyl phenyl sulfide

```json
{
    "name": "Or24a_methyl-phenyl-sulfide_complex",
    "modelSeeds": [1],
    "sequences": [
      {
        "protein": {
          "id": "A",
          "sequence": "precomputed MSA"
        }
      },
      {
        "ligand": {
          "id": "L",
          "smiles": "CSC1=CC=CC=C1"
        }
      }
    ],
    "dialect": "alphafold3",
    "version": 1
  }
```

<img src="pic/or24a-mps.png" alt="result" style="zoom:50%;" />

ipTM = 0.8, pTM = 0.71

### or85c and 2-heptanone

```json
{
    "name": "Or85c_2-heptanone_complex",
    "modelSeeds": [1],
    "sequences": [
      {
        "protein": {
          "id": "A",
          "sequence": "precomputed MSA"
        }
      },
      {
        "ligand": {
          "id": "L",
          "smiles": "CCCCCC(=O)C"
        }
      }
    ],
    "dialect": "alphafold3",
    "version": 1
  }
```

<img src="pic/or85c-2h.png" alt="result" style="zoom:50%;" />

ipTM = 0.94, pTM = 0.85

### or85c and 3-octanol (R)

```json
{
    "name": "Or85c_3-octanol-R_complex",
    "modelSeeds": [1],
    "sequences": [
      {
        "protein": {
          "id": "A",
          "sequence": "precomputed MSA"
        }
      },
      {
        "ligand": {
          "id": "L",
          "smiles": "CCCCCC[C@H](O)CC"
        }
      }
    ],
    "dialect": "alphafold3",
    "version": 1
  }
```

<img src="pic/or85c-3r.png" alt="result" style="zoom:50%;" />

ipTM = 0.94, pTM = 0.86

### or85c and 3-octanol (S)

```json
{
    "name": "Or85c_3-octanol-S_complex",
    "modelSeeds": [1],
    "sequences": [
      {
        "protein": {
          "id": "A",
          "sequence": "precomputed MSA"
        }
      },
      {
        "ligand": {
          "id": "L",
          "smiles": "CCCCCC(CC)O"
        }
      }
    ],
    "dialect": "alphafold3",
    "version": 1
  }
```

<img src="pic/or85c-3s.png" alt="result" style="zoom:50%;" />

ipTM = 0.94, pTM = 0.87