#!/usr/bin/bash

#SBATCH -J train-run
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-gpu=8
#SBATCH --mem-per-gpu=29G
#SBATCH -p batch_ce_ugrad
#SBATCH -t 1-0
#SBATCH -o /data/dydtkddhkdwk/20251021_MACE_FIT/logs/slurm-%A.out
#SBATCH -e /data/dydtkddhkdwk/20251021_MACE_FIT/logs/slurm-%A.err

# =======================
# 1. 환경 설정
# =======================
python CUDACHECK.py
echo "=== Job Completed ==="
nvidia-smi
date

