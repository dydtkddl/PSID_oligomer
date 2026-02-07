#!/usr/bin/bash

#SBATCH -J train-run
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-gpu=8
#SBATCH --mem-per-gpu=29G
#SBATCH -p batch_ce_ugrad
#SBATCH -t 1-0
#SBATCH -o /data/dydtkddhkdwk/20251021_MACE_FIT/logs/slurm-%A.out

# =======================
# 1. 환경 설정
# =======================
set -e  # 에러 시 종료
export CUDA_VISIBLE_DEVICES=0
export PYTHONUNBUFFERED=1  # 실시간 출력

# 작업 디렉토리 이동
cd /data/dydtkddhkdwk/20251021_MACE_FIT
mkdir -p logs

echo "=== MACE Training Started ==="
date
echo "Working Directory: $(pwd)"
nvidia-smi | head -n 20

# =======================
# 2. 실행
# =======================
python /data/dydtkddhkdwk/20251021_MACE_FIT/train_mace.py \
    --config /data/dydtkddhkdwk/20251021_MACE_FIT/config-02.yml \
    2>&1 | tee /data/dydtkddhkdwk/20251021_MACE_FIT/logs/train_$SLURM_JOB_ID.log

echo "=== Job Completed ==="
date
