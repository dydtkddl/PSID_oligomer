#!/usr/bin/bash

#SBATCH -J train-run
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-gpu=8
#SBATCH --mem-per-gpu=29G
#SBATCH -p batch_ce_ugrad
#SBATCH -t 1-0
#SBATCH -o /data/dydtkddhkdwk/20251023_MACE_MD_OPENMM/20251023_1305_BASICTEST/logs/slurm-%A.out
#SBATCH -e /data/dydtkddhkdwk/20251023_MACE_MD_OPENMM/20251023_1305_BASICTEST/logs/slurm-%A.err

# =======================
# 1. 환경 설정
# =======================
set -e  # 에러 시 종료
export CUDA_VISIBLE_DEVICES=0
export PYTHONUNBUFFERED=1  # 실시간 출력
export WORKDIR=$(pwd)
# 현재 작업 디렉토리를 기준으로 설정
WORKDIR=$(pwd)
cd "$WORKDIR"
mkdir -p logs

echo "=== MACE Training Started ==="
date
echo "Working Directory: $WORKDIR"
nvidia-smi | head -n 20

# =======================
# 2. 실행
# =======================
#sh "$WORKDIR/command_gpu.sh" \
#mace-md -f ejm_31.sdf --model_path MACE_test.model --output_dir md_test_2 --steps 1000 --unwrap --minimiser openmm \
mace-md -f ejm_31.sdf --model_path /data/dydtkddhkdwk/20251023_MACE_MD_OPENMM/20251023_1305_BASICTEST/MYMODEL/mace_off_finetune_stagetwo_compiled.model --output_dir md_test_2 --steps 1000 --unwrap --minimiser openmm \
    2>&1 | tee "$WORKDIR/logs/OPENMM_$SLURM_JOB_ID.log"

echo "=== Job Completed ==="
date

