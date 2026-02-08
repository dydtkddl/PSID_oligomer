#!/bin/bash
# install_cuda_13_1_wsl.sh

echo "=== CUDA Toolkit 13.1 WSL2 설치 ==="

# 기존 CUDA 정리
sudo apt-get --purge remove "*cublas*" "*cufft*" "*curand*" "*cusolver*" "*cusparse*" "*npp*" "*nvjpeg*" "cuda*" "nsight*" 2>/dev/null || true

# WSL-Ubuntu CUDA 저장소 추가
wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update

# CUDA Toolkit 13.1 설치 (최신 버전 확인 필요)
sudo apt-get install -y cuda-toolkit-12-3  # 13.1이 없으면 12.3 사용

# 환경 변수 설정
cat >> ~/.bashrc << 'EOF'

# CUDA 13.1 환경 설정 (RTX 5090)
export CUDA_HOME=/usr/local/cuda-12.3
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export CUDA_VISIBLE_DEVICES=0
EOF

source ~/.bashrc

# 설치 확인
echo ""
echo "=== CUDA 설치 확인 ==="
nvcc --version
nvidia-smi

echo "✅ CUDA Toolkit 설치 완료!"
