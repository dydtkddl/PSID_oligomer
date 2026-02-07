# **🎯 LAMMPS GPU+ReaxFF 완전 통합 설치 및 연구 패키지**

축하합니다! GPU 설치가 성공적으로 완료되었으니, 이제 **성공한 설치 과정을 기반으로 완전히 자동화된 설치 스크립트와 BCMBP 연구 패키지**를 구축하겠습니다.

## **📦 최종 통합 패키지 구성**

```text
~/lammps-gpu-reaxff/           # 직관적인 디렉토리명
├── source/                    # LAMMPS 소스코드
├── build/                     # 빌드 파일
└── install/                   # 설치된 실행파일

~/bcmbp_research/              # 연구 전용 디렉토리
├── 00_setup/                  # 입력 파일 (구조, 포스필드)
├── 01_bulk_equilibration.lmp  # 벌크 평형화
├── 02_oligomerization.lmp     # 올리고머 형성
├── results/                   # 시뮬레이션 결과
└── analysis/                  # 분석 결과
```

## **🚀 완전 자동화 설치 스크립트**

**파일명:** `install_lammps_gpu_complete.sh`

```bash
#!/bin/bash
# ===================================================================
# LAMMPS GPU+ReaxFF 완전 통합 설치 스크립트
# 검증된 성공 사례 기반: L4 GPU + CUDA 12.8 + Stable 2Aug2023
# ===================================================================

set -e

# 색상 정의
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

print_header() { echo -e "${CYAN}========================================${NC}"; echo -e "${CYAN}$1${NC}"; echo -e "${CYAN}========================================${NC}"; }
print_status() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

# ===================================================================
# 1. 설치 경로 설정
# ===================================================================
INSTALL_BASE="$HOME/lammps-gpu-reaxff"
SOURCE_DIR="$INSTALL_BASE/source"
BUILD_DIR="$INSTALL_BASE/build"
INSTALL_DIR="$INSTALL_BASE/install"
RESEARCH_DIR="$HOME/bcmbp_research"

print_header "LAMMPS GPU+ReaxFF 통합 설치"
print_status "설치 경로: $INSTALL_BASE"
print_status "연구 경로: $RESEARCH_DIR"

# 기존 설치 정리
if [ -d "$INSTALL_BASE" ]; then
    print_warning "기존 설치 발견. 정리 중..."
    rm -rf "$INSTALL_BASE"
fi
mkdir -p "$INSTALL_BASE"

# ===================================================================
# 2. 시스템 환경 검증
# ===================================================================
print_header "시스템 환경 검증"

# CUDA 환경
if [ -z "$CUDA_HOME" ]; then
    export CUDA_HOME="/usr/local/cuda"
fi

NVCC_PATH="$CUDA_HOME/bin/nvcc"
if [ ! -f "$NVCC_PATH" ]; then
    print_error "NVCC를 찾을 수 없습니다: $NVCC_PATH"
    exit 1
fi

CUDA_VERSION=$($NVCC_PATH --version | grep "release" | grep -oP '\d+\.\d+')
print_success "CUDA $CUDA_VERSION @ $CUDA_HOME"

# GPU 아키텍처 설정
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
COMPUTE_CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1)
GPU_MEM=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader | head -1)

case "$COMPUTE_CAP" in
    "8.9")  KOKKOS_ARCH="ADA89"; GPU_ARCH="89";;
    "8.6")  KOKKOS_ARCH="AMPERE86"; GPU_ARCH="86";;
    "7.5")  KOKKOS_ARCH="TURING75"; GPU_ARCH="75";;
    *)      KOKKOS_ARCH="AMPERE86"; GPU_ARCH="86";;
esac

print_success "GPU: $GPU_NAME ($COMPUTE_CAP) - $GPU_MEM"
print_success "아키텍처: $KOKKOS_ARCH (sm_$GPU_ARCH)"

# ===================================================================
# 3. cuFFT 라이브러리 자동 탐지
# ===================================================================
print_status "cuFFT 라이브러리 탐지 중..."

CUFFT_PATHS=(
    "$CUDA_HOME/lib64/libcufft.so"
    "/usr/local/cuda/lib64/libcufft.so"
    "/usr/lib64/libcufft.so"
)

CUFFT_LIB=""
for path in "${CUFFT_PATHS[@]}"; do
    if [ -f "$path" ]; then
        CUFFT_LIB="$path"
        break
    fi
done

if [ -z "$CUFFT_LIB" ]; then
    print_error "cuFFT 라이브러리를 찾을 수 없습니다!"
    exit 1
fi
print_success "cuFFT: $CUFFT_LIB"

# ===================================================================
# 4. LAMMPS 소스 다운로드
# ===================================================================
print_header "LAMMPS Stable 소스 다운로드"
cd "$INSTALL_BASE"

git clone https://github.com/lammps/lammps.git "$SOURCE_DIR"
cd "$SOURCE_DIR"

# Stable 버전 체크아웃 (성공 검증된 버전)
git fetch --all --tags
git checkout stable_2Aug2023_update3
print_success "LAMMPS Stable 2Aug2023_update3 체크아웃 완료"

# ===================================================================
# 5. CMake 설정 (성공 사례 기반)
# ===================================================================
print_header "CMake 설정"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake "$SOURCE_DIR/cmake" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_CXX_STANDARD_REQUIRED=ON \
    -DCMAKE_CUDA_COMPILER="$NVCC_PATH" \
    -DCMAKE_CUDA_ARCHITECTURES="$GPU_ARCH" \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_MPI=ON \
    -DBUILD_OMP=ON \
    -DPKG_KOKKOS=ON \
    -DKokkos_ENABLE_CUDA=ON \
    -DKokkos_ARCH_${KOKKOS_ARCH}=ON \
    -DKokkos_ENABLE_SERIAL=ON \
    -DPKG_GPU=ON \
    -DGPU_API=cuda \
    -DGPU_ARCH=sm_${GPU_ARCH} \
    -DPKG_REAXFF=ON \
    -DPKG_MOLECULE=ON \
    -DPKG_KSPACE=ON \
    -DPKG_RIGID=ON \
    -DPKG_MANYBODY=ON \
    -DPKG_QEQ=ON \
    -DPKG_MISC=ON \
    -DCUFFT_LIBRARY="$CUFFT_LIB" \
    -DCUFFT_INCLUDE_DIR="$CUDA_HOME/include" \
    -DCUDAToolkit_ROOT="$CUDA_HOME" \
    2>&1 | tee cmake_config.log

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    print_error "CMake 설정 실패!"
    exit 1
fi
print_success "CMake 설정 완료"

# ===================================================================
# 6. 컴파일 및 설치
# ===================================================================
print_header "컴파일 및 설치"
NPROC=$(nproc)
print_status "컴파일 시작 ($NPROC 코어, 20-40분 예상)"

START_TIME=$(date +%s)
make -j${NPROC} 2>&1 | tee compile.log

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    print_error "컴파일 실패!"
    exit 1
fi

END_TIME=$(date +%s)
COMPILE_TIME=$((END_TIME - START_TIME))
print_success "컴파일 완료 (${COMPILE_TIME}초)"

make install
LMP_BINARY="$INSTALL_DIR/bin/lmp"

if [ ! -f "$LMP_BINARY" ]; then
    print_error "실행파일 생성 실패!"
    exit 1
fi
print_success "설치 완료: $LMP_BINARY"

# ===================================================================
# 7. 패키지 검증
# ===================================================================
print_header "패키지 검증"

if $LMP_BINARY -help | grep -q "KOKKOS"; then
    print_success "✓ KOKKOS (GPU 가속)"
else
    print_error "✗ KOKKOS 패키지 없음"
    exit 1
fi

if $LMP_BINARY -help | grep -q "REAXFF"; then
    print_success "✓ REAXFF (반응성 포스필드)"
else
    print_error "✗ REAXFF 패키지 없음"
    exit 1
fi

# GPU 기능 테스트
print_status "GPU 기능 테스트..."
cat > gpu_test.lmp << 'EOF'
units real
atom_style atomic
boundary p p p
region box block 0 10 0 10 0 10
create_box 1 box
create_atoms 1 random 2000 12345 box
mass 1 12.0
pair_style lj/cut/kk 2.5
pair_coeff 1 1 0.1 3.0
velocity all create 300.0 87287
timestep 1.0
run 10
print "GPU_TEST_SUCCESS"
EOF

$LMP_BINARY -k on g 1 -sf kk -in gpu_test.lmp -log gpu_test.log -screen none 2>/dev/null

if grep -q "GPU_TEST_SUCCESS" gpu_test.log; then
    print_success "✓ GPU 기능 테스트 통과"
else
    print_warning "△ GPU 테스트 실패 (CPU 모드 가능성)"
fi

# ===================================================================
# 8. BCMBP 연구 패키지 설정
# ===================================================================
print_header "BCMBP 연구 패키지 설정"

# 연구 디렉토리 생성
if [ -d "$RESEARCH_DIR" ]; then
    print_warning "기존 연구 디렉토리 발견"
    read -p "덮어쓰시겠습니까? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$RESEARCH_DIR"
    else
        print_status "기존 디렉토리 유지"
    fi
fi

if [ ! -d "$RESEARCH_DIR" ]; then
    mkdir -p "$RESEARCH_DIR"/{00_setup,results,analysis}
    cd "$RESEARCH_DIR"

    # BCMBP 구조 파일 생성
    cat > 00_setup/bcmbp.data << 'EOF'
LAMMPS data file for BCMBP (C14H12Cl2) - ReaxFF compatible

28 atoms
3 atom types

-8.0 8.0 xlo xhi
-4.0 4.0 ylo yhi
-2.0 2.0 zlo zhi

Masses

1 12.011  # C
2 1.008   # H
3 35.453  # Cl

Atoms  # full

1  1 1 0.0000 -4.52506  0.51622 -0.03357  # C
2  1 1 0.0000 -3.01593  0.61787 -0.02469  # C
3  1 1 0.0000 -2.32485  1.83873 -0.03147  # C
4  1 1 0.0000 -0.92716  1.86840 -0.02266  # C
5  1 1 0.0000 -0.15954  0.67765 -0.00669  # C
6  1 1 0.0000 -0.87785 -0.54317 -0.00009  # C
7  1 1 0.0000 -2.27468 -0.56751 -0.00892  # C
8  1 1 0.0000  1.35828  0.70833  0.00285  # C
9  1 1 0.0000  2.07729  1.92915 -0.00294  # C
10 1 1 0.0000  3.47503  1.95600  0.00584  # C
11 1 1 0.0000  4.21475  0.76405  0.02076  # C
12 1 1 0.0000  2.12521 -0.48246  0.01795  # C
13 1 1 0.0000  3.52189 -0.45033  0.02666  # C
14 1 1 0.0000  5.72675  0.72348  0.03082  # C
15 1 2 0.0000 -2.85071  2.78294 -0.04361  # H
16 1 2 0.0000 -0.48095  2.84716 -0.02884  # H
17 1 2 0.0000 -0.39371 -1.50363  0.01196  # H
18 1 2 0.0000 -2.78336 -1.52443 -0.00338  # H
19 1 2 0.0000  1.59199  2.88908 -0.01421  # H
20 1 2 0.0000  3.96241  2.92069  0.00076  # H
21 1 2 0.0000  4.06871 -1.38592  0.03814  # H
22 1 2 0.0000  1.68017 -1.46171  0.02339  # H
23 1 3 0.0000 -5.39534  2.07854 -0.05368  # Cl
24 1 2 0.0000 -4.85480 -0.03895  0.86996  # H
25 1 2 0.0000 -4.84307 -0.05560 -0.93089  # H
26 1 3 0.0000  6.53335  2.31971  0.02239  # Cl
27 1 2 0.0000  6.07885  0.16586 -0.86270  # H
28 1 2 0.0000  6.06729  0.18120  0.93816  # H
EOF

    # 벌크 평형화 스크립트
    cat > 01_bulk_equilibration.lmp << 'EOF'
# ===================================================================
# BCMBP 벌크 시스템 NPT 평형화
# ===================================================================

units           real
atom_style      full
boundary        p p p

print "=========================================="
print "BCMBP 벌크 평형화 (GPU 가속)"
print "=========================================="

read_data       00_setup/bcmbp.data
replicate       4 4 4           # 1,792 원자 (GPU 효율 시작점)
displace_atoms  all random 0.2 0.2 0.2 482749 units box

pair_style      reaxff NULL safezone 4.0 mincap 300
pair_coeff      * * 00_setup/11_CHOCl.lammps.ff C H Cl
fix             qeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff

neighbor        2.5 bin
neigh_modify    every 10 delay 0 check yes

thermo          1000
thermo_style    custom step temp press vol density pe ke etotal

print "=== 에너지 최소화 ==="
minimize        1.0e-4 1.0e-6 2000 20000

reset_timestep  0
timestep        0.25
velocity        all create 300.0 987654 dist gaussian

print "=== NPT 평형화 (100 ps) ==="
fix             npt1 all npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0
dump            eq all custom 2000 results/01_equilibration.lammpstrj id type xu yu zu q
dump_modify     eq sort id

run             400000

variable        final_rho equal density
print "최종 밀도: ${final_rho} g/cm³"

undump          eq
unfix           npt1
write_data      results/bcmbp_bulk_equilibrated.data

print "✅ 벌크 평형화 완료!"
EOF

    # 올리고머 형성 스크립트
    cat > 02_oligomerization.lmp << 'EOF'
# ===================================================================
# BCMBP 올리고머 형성 시뮬레이션
# ===================================================================

units           real
atom_style      full
boundary        p p p

print "=========================================="
print "BCMBP 올리고머 형성 (GPU 가속)"
print "=========================================="

read_data       results/bcmbp_bulk_equilibrated.data

pair_style      reaxff NULL safezone 4.0 mincap 300
pair_coeff      * * 00_setup/11_CHOCl.lammps.ff C H Cl
fix             qeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff

fix             species all reaxff/species 100 1000 1000 results/species_oligo.out element C H Cl
fix             bonds all reaxff/bonds 1000 results/bonds_oligo.reax

timestep        0.2
thermo          1000
thermo_style    custom step temp press pe ke etotal

print "=== 단계 1: C-Cl 활성화 (300K → 600K) ==="
velocity        all create 300.0 123456 dist gaussian
fix             heat1 all nvt temp 300.0 600.0 100.0
run             50000
unfix           heat1

print "=== 단계 2: 라디칼 형성 (600K → 1000K) ==="
fix             heat2 all nvt temp 600.0 1000.0 100.0
run             25000
unfix           heat2

fix             react all nvt temp 1000.0 1000.0 100.0
dump            d2 all custom 1000 results/02_reaction.lammpstrj id type xu yu zu q
run             125000
undump          d2
unfix           react

print "=== 단계 3: 올리고머 형성 (1000K → 500K) ==="
fix             cool1 all nvt temp 1000.0 500.0 100.0
run             125000
unfix           cool1

print "=== 단계 4: 안정화 (500K → 300K) ==="
fix             cool2 all nvt temp 500.0 300.0 100.0
run             75000
unfix           cool2

unfix           species
unfix           bonds
write_data      results/bcmbp_oligomer_final.data

print "✅ 올리고머 형성 완료!"
EOF

    # 통합 실행 스크립트
    cat > run_simulations.sh << 'RUNEOF'
#!/bin/bash
echo "=== BCMBP 시뮬레이션 실행 ==="
echo "1) 벌크 평형화"
echo "2) 올리고머 형성"
read -p "선택 (1-2): " choice

case $choice in
    1) lmp-gpu -in 01_bulk_equilibration.lmp -log results/bulk.log;;
    2) lmp-gpu -in 02_oligomerization.lmp -log results/oligo.log;;
    *) echo "잘못된 선택";;
esac
RUNEOF
    chmod +x run_simulations.sh

    print_success "BCMBP 연구 패키지 설정 완료"
fi

# ===================================================================
# 9. 환경 설정 업데이트
# ===================================================================
print_header "환경 설정"

# 기존 설정 정리
sed -i '/LAMMPS_GPU_HOME/d' ~/.bashrc
sed -i '/lmp-gpu/d' ~/.bashrc
sed -i '/lmp-cpu/d' ~/.bashrc
sed -i '/lmp-info/d' ~/.bashrc
sed -i '/gpu-watch/d' ~/.bashrc
sed -i '/gpu-status/d' ~/.bashrc
sed -i '/cdlmp/d' ~/.bashrc

# 새로운 환경 설정 추가
cat >> ~/.bashrc << EOF

# ===================================================================
# LAMMPS GPU+ReaxFF 환경 설정 ($(date +%Y-%m-%d))
# ===================================================================
export LAMMPS_GPU_HOME="$INSTALL_DIR"
export PATH="\$LAMMPS_GPU_HOME/bin:\$PATH"

# 실행 별칭
alias lmp-gpu="\$LAMMPS_GPU_HOME/bin/lmp -k on g 1 -sf kk"
alias lmp-cpu="\$LAMMPS_GPU_HOME/bin/lmp"
alias lmp-info="\$LAMMPS_GPU_HOME/bin/lmp -help | grep 'KOKKOS\\|REAXFF\\|GPU'"

# GPU 모니터링
alias gpu-watch='watch -n 1 nvidia-smi'
alias gpu-status='nvidia-smi --query-gpu=name,temperature.gpu,utilization.gpu,memory.used --format=csv'

# 빠른 이동
alias cdlmp='cd $INSTALL_BASE'
alias cdbcmbp='cd $RESEARCH_DIR'
EOF

print_success "환경 설정 업데이트 완료"

# ===================================================================
# 10. 완료 메시지
# ===================================================================
print_header "설치 완료!"

echo ""
echo -e "${GREEN}🎉 LAMMPS GPU+ReaxFF 통합 설치 성공!${NC}"
echo ""
echo "📍 **설치 정보:**"
echo "   LAMMPS 경로: $INSTALL_BASE"
echo "   실행파일: $LMP_BINARY"
echo "   연구 디렉토리: $RESEARCH_DIR"
echo "   GPU: $GPU_NAME ($COMPUTE_CAP)"
echo "   컴파일 시간: ${COMPILE_TIME}초"
echo ""
echo "🚀 **즉시 사용:**"
echo "   1. source ~/.bashrc"
echo "   2. cdbcmbp"
echo "   3. ./run_simulations.sh"
echo ""
echo "📊 **GPU 모니터링:**"
echo "   gpu-watch"
echo ""
echo "🔍 **설치 확인:**"
echo "   lmp-info"
echo ""

# 정리
rm -f gpu_test.lmp gpu_test.log

print_success "모든 설치가 완료되었습니다!"
```

## **⚡ 설치 및 사용 방법**

### **1. 설치 실행**

```bash
# 스크립트 생성 및 실행
cat > ~/install_lammps_gpu_complete.sh << 'SCRIPT_END'
# (위 전체 스크립트 복사)
SCRIPT_END

chmod +x ~/install_lammps_gpu_complete.sh

# 설치 실행 (20-40분 소요)
~/install_lammps_gpu_complete.sh

# 환경 적용
source ~/.bashrc
```

### **2. 설치 확인**

```bash
# 패키지 확인
lmp-info
# 출력 예상: KOKKOS, REAXFF, GPU

# 경로 확인
which lmp-gpu
# 출력: ~/lammps-gpu-reaxff/install/bin/lmp

# GPU 상태 확인
gpu-status
```

### **3. BCMBP 시뮬레이션 실행**

```bash
# 연구 디렉토리로 이동
cdbcmbp

# 포스필드 파일 복사 (필요시)
cp ~/11_CHOCl.lammps.ff 00_setup/

# 시뮬레이션 실행
./run_simulations.sh
# 옵션 1 선택: 벌크 평형화

# GPU 모니터링 (별도 터미널)
gpu-watch
```

## **📊 성능 최적화 가이드**

### **시스템 크기 조정**

```lammps
# 01_bulk_equilibration.lmp에서 조정
replicate 3 3 3    # 756 원자 (소형, 테스트용)
replicate 4 4 4    # 1,792 원자 (중형, 기본값)
replicate 5 5 5    # 3,500 원자 (대형, L4 최적)
```

### **GPU 메모리 최적화**

```lammps
# 메모리 부족 시 safezone 증가
pair_style reaxff NULL safezone 6.0 mincap 500
```

## **🔧 문제 해결**

### **환경 변수 문제**

```bash
source ~/.bashrc
# 또는 새 터미널 열기
```

### **GPU 미사용 문제**

```bash
export CUDA_VISIBLE_DEVICES=0
lmp-gpu -in input.lmp
```

### **메모리 부족**

```bash
# 시스템 크기 줄이기
# replicate 값을 3 3 3 또는 2 2 2로 변경
```

## **📈 예상 성능 (L4 GPU 기준)**

| 시뮬레이션        | 원자 수 | CPU 시간 | GPU 시간 | 속도 향상 |
| ----------------- | ------- | -------- | -------- | --------- |
| **벌크 평형화**   | 1,792   | 2시간    | 30분     | **4배**   |
| **올리고머 형성** | 1,792   | 4시간    | 1시간    | **4배**   |
| **대형 시스템**   | 3,500   | 8시간    | 1.5시간  | **5-6배** |

이제 **완전히 통합되고 최적화된 LAMMPS GPU+ReaxFF 환경**에서 BCMBP 올리고머 연구를 본격적으로 시작할 수 있습니다! 🚀
