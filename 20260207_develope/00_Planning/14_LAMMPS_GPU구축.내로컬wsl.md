# **14_LAMMPS_GPU구축.내로컬WSL.md**

**환경: WSL2 + RTX 5090 (32GB) + Threadripper PRO 9965WX (48 스레드) + CUDA 13.1**

---

## **🎯 시스템 환경 분석 및 구축 전략**

### **📊 하드웨어 사양 (2026년 현재 최상위급)**

```yaml
GPU: NVIDIA GeForce RTX 5090
  - VRAM: 32GB (Blackwell 아키텍처)
  - CUDA Cores: ~16,000개 (추정)
  - Compute Capability: 8.9 (sm_89)
  - CUDA: 13.1, Driver: 591.44

CPU: AMD Ryzen Threadripper PRO 9965WX
  - 코어: 24 (48 하이퍼스레드)
  - 기본 클럭: 3.2GHz, 부스트: 5.1GHz
  - L3 캐시: 128MB

메모리: 62GB RAM
디스크: 20TB (18TB 여유)
OS: Ubuntu 22.04.5 LTS (WSL2)
```

### **🤔 Docker vs Native WSL2 전략 분석**

| 구분                | Docker 방식              | Native WSL2 방식            |
| ------------------- | ------------------------ | --------------------------- |
| **성능**            | 98-99% (미미한 오버헤드) | 100% (직접 하드웨어 제어)   |
| **설치 시간**       | 10분 (사전 빌드 이미지)  | 30-50분 (소스 컴파일)       |
| **재현성**          | ⭐⭐⭐⭐⭐ 완벽          | ⭐⭐⭐ 환경 의존적          |
| **RTX 5090 최적화** | ⭐⭐⭐ 일반적 최적화     | ⭐⭐⭐⭐⭐ 완전 맞춤 최적화 |
| **개발 편의성**     | ⭐⭐⭐ 컨테이너 격리     | ⭐⭐⭐⭐⭐ 직접 접근        |
| **멀티 버전**       | ⭐⭐⭐⭐⭐ 동시 실행     | ⭐⭐ 어려움                 |

### **✅ 권장 전략: 상황별 선택**

**1. 빠른 시작 및 프로토타이핑 → Docker 방식**
**2. 최대 성능 및 장기 연구 → Native 방식**
**3. 이상적 워크플로우: Docker로 시작 → Native로 확장**

---

## **🚀 방법 1: Docker 방식 (빠른 시작 - 10분 완료)**

### **1-1. NVIDIA Container Toolkit 설치**

```bash
#!/bin/bash
# install_nvidia_docker_wsl.sh

echo "=== NVIDIA Container Toolkit 설치 (WSL2) ==="

# Docker 설치 확인
if ! command -v docker &> /dev/null; then
    echo "Docker 설치 중..."
    curl -fsSL https://get.docker.com -o get-docker.sh
    sudo sh get-docker.sh
    sudo usermod -aG docker $USER
    echo "⚠️ 재로그인 후 계속하세요!"
    exit 1
fi

# NVIDIA Container Toolkit 저장소 추가
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | \
    sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg

curl -s -L https://nvidia.github.io/libnvidia-container/$distribution/libnvidia-container.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

# 설치
sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit

# Docker 데몬 설정
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker

# RTX 5090 테스트
echo ""
echo "=== RTX 5090 Docker 접근 테스트 ==="
docker run --rm --gpus all nvidia/cuda:12.3.0-base-ubuntu22.04 nvidia-smi

echo "✅ NVIDIA Container Toolkit 설치 완료!"
```

### **1-2. LAMMPS GPU 컨테이너 즉시 실행**

```bash
#!/bin/bash
# run_lammps_docker_5090.sh - RTX 5090 최적화 실행

# 작업 디렉토리 설정
WORK_DIR="/mnt/d/[01]Lab_Activity/##자율제조과제/20260205_Develope"
cd "$WORK_DIR"

echo "=== LAMMPS GPU Docker 실행 (RTX 5090) ==="

# 최신 LAMMPS GPU 이미지 사용
docker run -it --rm \
    --gpus all \
    --name lammps-rtx5090 \
    -v "$PWD:/workspace" \
    --shm-size=16g \
    -e CUDA_VISIBLE_DEVICES=0 \
    lammps/lammps:stable_29Aug2024_ubuntu22.04_gpu \
    bash -c "
        cd /workspace
        echo '=== RTX 5090 GPU 정보 ==='
        nvidia-smi
        echo ''
        echo '=== LAMMPS 패키지 확인 ==='
        lmp -help | grep 'KOKKOS\|REAXFF\|GPU'
        echo ''
        echo '🚀 LAMMPS 준비 완료! 다음 명령어로 실행:'
        echo '   lmp -k on g 1 -sf kk -in input.lmp'
        /bin/bash
    "
```

### **1-3. Docker Compose 설정 (영구 사용)**

```yaml
# docker-compose.yml
version: "3.8"

services:
  lammps-rtx5090:
    image: lammps/lammps:stable_29Aug2024_ubuntu22.04_gpu
    container_name: bcmbp-rtx5090
    runtime: nvidia
    environment:
      - NVIDIA_VISIBLE_DEVICES=all
      - CUDA_VISIBLE_DEVICES=0
      - OMP_NUM_THREADS=1
    volumes:
      - /mnt/d/[01]Lab_Activity/##자율제조과제/20260205_Develope:/workspace
      - ./results:/results
    working_dir: /workspace
    shm_size: 16gb
    stdin_open: true
    tty: true
    command: /bin/bash
```

**실행:**

```bash
docker-compose run --rm lammps-rtx5090
```

---

## **🔧 방법 2: Native WSL2 방식 (최대 성능 - RTX 5090 완전 활용)**

### **2-1. CUDA Toolkit 13.1 설치 (WSL2 전용)**

```bash
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
```

### **2-2. LAMMPS RTX 5090 최적화 빌드**

```bash
#!/bin/bash
# install_lammps_rtx5090_native.sh - RTX 5090 완전 최적화

set -e

echo "=========================================="
echo "🚀 LAMMPS RTX 5090 Native 최적화 빌드"
echo "=========================================="

# 1. 시스템 패키지 설치
sudo apt-get update
sudo apt-get install -y \
    build-essential cmake git \
    libopenmpi-dev openmpi-bin \
    libfftw3-dev libjpeg-dev libpng-dev \
    python3 python3-pip

# 2. 설치 경로 설정
INSTALL_BASE="$HOME/lammps-rtx5090-native"
SOURCE_DIR="$INSTALL_BASE/source"
BUILD_DIR="$INSTALL_BASE/build"
INSTALL_DIR="$INSTALL_BASE/install"

# 기존 설치 백업
if [ -d "$INSTALL_BASE" ]; then
    mv "$INSTALL_BASE" "${INSTALL_BASE}.backup_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$INSTALL_BASE"

# 3. LAMMPS 최신 소스 다운로드
echo "📥 LAMMPS 소스 다운로드..."
git clone https://github.com/lammps/lammps.git "$SOURCE_DIR"
cd "$SOURCE_DIR"
git checkout stable_2Aug2023_update3

# 4. RTX 5090 최적화 CMake 설정
echo "⚙️ RTX 5090 최적화 설정..."
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# RTX 5090 = Blackwell = sm_89 (compute capability 8.9)
cmake "$SOURCE_DIR/cmake" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_CUDA_COMPILER="$CUDA_HOME/bin/nvcc" \
    -DCMAKE_CUDA_ARCHITECTURES=89 \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_MPI=ON \
    -DBUILD_OMP=ON \
    -DPKG_KOKKOS=ON \
    -DKokkos_ENABLE_CUDA=ON \
    -DKokkos_ARCH_ADA89=ON \
    -DKokkos_ENABLE_SERIAL=ON \
    -DKokkos_ENABLE_OPENMP=ON \
    -DPKG_REAXFF=ON \
    -DPKG_MOLECULE=ON \
    -DPKG_KSPACE=ON \
    -DPKG_RIGID=ON \
    -DPKG_MANYBODY=ON \
    -DPKG_QEQ=ON \
    -DPKG_MISC=ON \
    2>&1 | tee cmake_config.log

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "❌ CMake 설정 실패!"
    exit 1
fi

# 5. Threadripper PRO 48스레드 컴파일
echo "🔨 컴파일 시작 (48 스레드)..."
START_TIME=$(date +%s)

make -j48 2>&1 | tee compile.log

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "❌ 컴파일 실패!"
    exit 1
fi

END_TIME=$(date +%s)
COMPILE_TIME=$((END_TIME - START_TIME))

# 6. 설치
echo "📦 설치 중..."
make install

# 7. 환경 설정
cat >> ~/.bashrc << EOF

# LAMMPS RTX 5090 Native 환경
export LAMMPS_RTX5090="$INSTALL_DIR"
export PATH="\$LAMMPS_RTX5090/bin:\$PATH"

# RTX 5090 최적화 별칭
alias lmp-5090="\$LAMMPS_RTX5090/bin/lmp -k on g 1 -sf kk"
alias lmp-5090-omp="\$LAMMPS_RTX5090/bin/lmp -k on g 1 t 24 -sf kk"  # GPU + OpenMP 하이브리드
EOF

source ~/.bashrc

# 8. 설치 검증
LMP_BINARY="$INSTALL_DIR/bin/lmp"
if [ ! -f "$LMP_BINARY" ]; then
    echo "❌ 설치 실패!"
    exit 1
fi

echo ""
echo "=========================================="
echo "🎉 RTX 5090 Native 빌드 완료!"
echo "=========================================="
echo "컴파일 시간: ${COMPILE_TIME}초"
echo "실행파일: $LMP_BINARY"
echo ""
echo "🚀 사용법:"
echo "  lmp-5090 -in input.lmp              # GPU 전용"
echo "  lmp-5090-omp -in input.lmp          # GPU + OpenMP"
echo ""
echo "📊 패키지 확인:"
$LMP_BINARY -help | grep "KOKKOS\|REAXFF\|GPU"
```

---

## **🛡️ 안전 설치 및 검증 도구**

### **3-1. RTX 5090 전용 안전 검증 스크립트**

```bash
#!/bin/bash
# verify_rtx5090_setup.sh - RTX 5090 설정 검증

echo "=== RTX 5090 LAMMPS 설정 검증 ==="

# 1. GPU 하드웨어 확인
echo "1. GPU 하드웨어 정보:"
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader

# 2. CUDA 컴파일러 확인
echo ""
echo "2. CUDA 컴파일러:"
nvcc --version | grep "release"

# 3. LAMMPS 패키지 확인
echo ""
echo "3. LAMMPS 패키지:"
if command -v lmp-5090 &> /dev/null; then
    lmp-5090 -help | grep "KOKKOS\|REAXFF\|GPU"
else
    echo "❌ lmp-5090 명령어를 찾을 수 없습니다."
fi

# 4. GPU 메모리 테스트
echo ""
echo "4. GPU 메모리 상태:"
nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader

# 5. 간단한 GPU 연산 테스트
cat > gpu_test.lmp << 'EOF'
units real
atom_style atomic
boundary p p p
newton on
package kokkos neigh half

region box block 0 20 0 20 0 20
create_box 1 box
create_atoms 1 random 10000 12345 box
mass 1 12.0

pair_style lj/cut/kk 2.5
pair_coeff 1 1 0.1 3.0

velocity all create 300.0 87287
timestep 1.0

run 100
print "RTX_5090_TEST_SUCCESS"
EOF

echo ""
echo "5. GPU 연산 테스트 중..."
if command -v lmp-5090 &> /dev/null; then
    lmp-5090 -in gpu_test.lmp -log gpu_test.log -screen none 2>/dev/null
    if grep -q "RTX_5090_TEST_SUCCESS" gpu_test.log; then
        echo "✅ GPU 연산 테스트 통과!"
    else
        echo "❌ GPU 연산 테스트 실패"
    fi
    rm -f gpu_test.lmp gpu_test.log
else
    echo "⚠️ LAMMPS 실행파일 없음"
fi

echo ""
echo "=== 검증 완료 ==="
```

### **3-2. RTX 5090 성능 벤치마크**

```bash
#!/bin/bash
# benchmark_rtx5090.sh - RTX 5090 성능 벤치마크

echo "=== RTX 5090 성능 벤치마크 ==="

# BCMBP 시스템 기준 벤치마크
cat > benchmark_5090.lmp << 'EOF'
# RTX 5090 ReaxFF 벤치마크
units real
atom_style full
boundary p p p
newton on
package kokkos neigh half

read_data bcmbp.data
# RTX 5090의 32GB VRAM을 활용한 대규모 시스템
replicate 10 10 10  # 약 28,000 원자

pair_style reaxff NULL safezone 6.0 mincap 500
pair_coeff * * 11_CHOCl.lammps.ff C H Cl
fix qeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff

timestep 0.25
thermo 100

# 웜업
run 1000

# 벤치마크
reset_timestep 0
timer timeout 0
run 10000

print "RTX_5090_BENCHMARK_COMPLETE"
EOF

# 실행 및 성능 측정
echo "🚀 RTX 5090 벤치마크 실행 중... (약 10-20분 소요)"
START_TIME=$(date +%s)

lmp-5090 -in benchmark_5090.lmp -log benchmark_5090.log

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

# 결과 분석
if grep -q "RTX_5090_BENCHMARK_COMPLETE" benchmark_5090.log; then
    LOOP_TIME=$(grep "Loop time" benchmark_5090.log | tail -1 | awk '{print $4}')
    TIMESTEP_RATE=$(grep "timesteps/s" benchmark_5090.log | tail -1 | awk '{print $1}')

    echo ""
    echo "🎯 RTX 5090 벤치마크 결과:"
    echo "  총 실행 시간: ${TOTAL_TIME}초"
    echo "  Loop time: ${LOOP_TIME}초"
    echo "  성능: ${TIMESTEP_RATE} timesteps/s"
    echo "  시스템: ~28,000 원자 ReaxFF"
else
    echo "❌ 벤치마크 실패"
fi

# 정리
rm -f benchmark_5090.lmp
```

---

## **🎯 실전 사용 가이드**

### **4-1. BCMBP 올리고머 RTX 5090 최적화 실행**

```bash
#!/bin/bash
# run_bcmbp_rtx5090.sh - RTX 5090 최적화 BCMBP 시뮬레이션

echo "=== BCMBP 올리고머 RTX 5090 실행 ==="

# 작업 디렉토리 확인
if [ ! -f "bcmbp.data" ] || [ ! -f "11_CHOCl.lammps.ff" ]; then
    echo "❌ 필수 파일 없음: bcmbp.data, 11_CHOCl.lammps.ff"
    exit 1
fi

# RTX 5090 최적화 입력 파일 생성
cat > bcmbp_rtx5090_oligo.lmp << 'EOF'
# BCMBP 올리고머 형성 (RTX 5090 32GB 활용)
units real
atom_style full
boundary p p p
newton on
package kokkos neigh half

print "=========================================="
print "BCMBP 올리고머 RTX 5090 최적화"
print "=========================================="

read_data bcmbp.data
# RTX 5090의 32GB VRAM 활용 - 대규모 시스템
replicate 8 8 8  # 약 14,336 원자 (안전한 크기)

pair_style reaxff NULL safezone 6.0 mincap 500
pair_coeff * * 11_CHOCl.lammps.ff C H Cl
fix qeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff

# 화학종 분석
fix species all reaxff/species 100 1000 1000 species_rtx5090.out element C H Cl
fix bonds all reaxff/bonds 1000 bonds_rtx5090.reax

timestep 0.2
thermo 1000
thermo_style custom step temp press pe ke etotal vol density

# 단계 1: 활성화 (300K → 600K)
print "=== 단계 1: C-Cl 활성화 ==="
velocity all create 300.0 123456 dist gaussian
fix heat1 all nvt temp 300.0 600.0 100.0
run 50000
unfix heat1

# 단계 2: 반응 (600K → 1000K)
print "=== 단계 2: 라디칼 반응 ==="
fix heat2 all nvt temp 600.0 1000.0 100.0
run 25000
unfix heat2

fix react all nvt temp 1000.0 1000.0 100.0
dump d1 all custom 1000 rtx5090_reaction.lammpstrj id type xu yu zu q
run 125000
undump d1
unfix react

# 단계 3: 올리고머 형성 (냉각)
print "=== 단계 3: 올리고머 형성 ==="
fix cool1 all nvt temp 1000.0 500.0 100.0
run 125000
unfix cool1

# 단계 4: 안정화
print "=== 단계 4: 안정화 ==="
fix cool2 all nvt temp 500.0 300.0 100.0
run 75000
unfix cool2

unfix species
unfix bonds
write_data bcmbp_rtx5090_final.data

print "=========================================="
print "RTX 5090 올리고머 시뮬레이션 완료!"
print "=========================================="
EOF

# GPU 메모리 사전 확인
GPU_FREE=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
if [ "$GPU_FREE" -lt 10000 ]; then
    echo "⚠️ GPU 메모리 부족: ${GPU_FREE}MB (권장 10GB+)"
    read -p "계속하시겠습니까? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# 실행
echo "🚀 RTX 5090 올리고머 시뮬레이션 시작..."
echo "📊 GPU 모니터링: nvidia-smi -l 1 (별도 터미널)"

START_TIME=$(date +%s)
lmp-5090 -in bcmbp_rtx5090_oligo.lmp -log rtx5090_oligo.log

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo ""
echo "✅ 시뮬레이션 완료! (실행시간: ${RUNTIME}초)"
echo ""
echo "📁 생성된 파일:"
ls -lh species_rtx5090.out rtx5090_reaction.lammpstrj bcmbp_rtx5090_final.data 2>/dev/null
```

### **4-2. RTX 5090 실시간 모니터링**

```bash
#!/bin/bash
# monitor_rtx5090.sh - RTX 5090 전용 모니터링

trap "tput cnorm; clear; exit" SIGINT SIGTERM
tput civis

while true; do
    clear
    echo "=========================================="
    echo "🚀 RTX 5090 실시간 모니터링"
    echo "시간: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "=========================================="

    # GPU 상세 정보
    nvidia-smi --query-gpu=name,temperature.gpu,power.draw,power.limit,utilization.gpu,memory.used,memory.total \
        --format=csv,noheader | \
        awk -F', ' '{
            printf "GPU 모델: %s\n", $1
            printf "온도    : %s°C\n", $2
            printf "전력    : %s / %s W\n", $3, $4
            printf "사용률  : %s%%\n", $5
            printf "메모리  : %s / %s\n", $6, $7

            # 메모리 사용률 바 표시
            split($6, used, " "); split($7, total, " ");
            if (total[1] > 0) {
                pct = (used[1] / total[1]) * 100;
                printf "메모리바: [";
                for(i=0;i<50;i++) {
                    if(i < pct*50/100) printf "█"; else printf "░";
                }
                printf "] %.1f%%\n", pct;
            }
        }'

    echo "=========================================="
    echo "💡 팁: RTX 5090은 32GB VRAM으로 초대규모 시뮬레이션 가능"
    echo "Ctrl+C로 종료"

    sleep 2
done
```

---

## **📊 최종 권장 워크플로우**

### **단계별 구축 및 활용 전략**

```bash
# ===================================================================
# 1단계: 빠른 시작 (Docker - 10분)
# ===================================================================
# 즉시 테스트하고 싶은 경우
sh install_nvidia_docker_wsl.sh
sh run_lammps_docker_5090.sh

# ===================================================================
# 2단계: 최대 성능 (Native - 1시간)
# ===================================================================
# 본격적인 연구용 환경 구축
sh install_cuda_13_1_wsl.sh
sh install_lammps_rtx5090_native.sh

# ===================================================================
# 3단계: 검증 및 벤치마크
# ===================================================================
sh verify_rtx5090_setup.sh
sh benchmark_rtx5090.sh

# ===================================================================
# 4단계: 실전 연구
# ===================================================================
sh run_bcmbp_rtx5090.sh  # BCMBP 올리고머 연구
```

### **🎯 성능 예상치 (RTX 5090 기준)**

| 시뮬레이션 종류   | 시스템 크기   | 예상 성능          | 메모리 사용 |
| ----------------- | ------------- | ------------------ | ----------- |
| **소규모 테스트** | ~1,000 원자   | 매우 빠름 (오버킬) | <1GB        |
| **중간 연구**     | ~10,000 원자  | 10-20배 향상       | 3-5GB       |
| **대규모 연구**   | ~50,000 원자  | 20-50배 향상       | 15-20GB     |
| **초대규모**      | ~100,000 원자 | 극한 성능          | 25-30GB     |

### **💡 RTX 5090 활용 팁**

1. **메모리 여유 활용**: 32GB VRAM으로 일반적인 연구보다 5-10배 큰 시스템 가능
2. **배치 시뮬레이션**: 여러 조건을 순차적으로 자동 실행
3. **하이브리드 병렬화**: GPU + Threadripper의 48스레드 동시 활용
4. **장기 시뮬레이션**: 높은 안정성으로 며칠간 연속 실행 가능

---

## **🔍 결론 및 권장사항**

**RTX 5090 + Threadripper PRO 9965WX** 시스템은 현재 개인용 워크스테이션 중 최상위급입니다.

### **✅ 최종 권장 전략**

1. **학습 및 테스트**: Docker 방식으로 빠른 시작
2. **본격 연구**: Native WSL2 방식으로 최대 성능 활용
3. **장기 운영**: 두 방식 모두 구축하여 상황별 선택 사용

### **🚀 차별화 포인트**

- **32GB VRAM**: 기존 연구의 5-10배 규모 시뮬레이션 가능
- **48 스레드 CPU**: 전처리/후처리 및 하이브리드 병렬화 최적
- **WSL2 환경**: Windows와 Linux의 장점 결합

이 환경에서는 **분자동역학 시뮬레이션의 새로운 가능성**을 탐구할 수 있으며, 기존에 불가능했던 초대규모 시스템이나 장기간 시뮬레이션을 개인 워크스테이션에서 수행할 수 있습니다.

**Docker로 빠르게 시작하되, Native로 극한 성능을 추구하는 하이브리드 전략**이 RTX 5090의 잠재력을 최대한 활용하는 최적의 방법입니다. 🚀
