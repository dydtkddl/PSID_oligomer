#!/bin/bash
export OMP_NUM_THREADS=30

PART=5  # ← 여기 숫자만 바꾸면 됨

python 02_xtb_evaluate.py \
    --input 01_all_xtb_sampled/01_all_xtb_sampled_part_${PART}.xyz \
    --out_xyz ionpairs_xtb_sample_${PART}.xyz \
    --out_npz ionpairs_xtb_sample_${PART}.npz \
    --method GFN2-xTB

