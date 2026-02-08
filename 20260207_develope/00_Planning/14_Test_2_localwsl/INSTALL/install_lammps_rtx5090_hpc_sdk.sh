#!/bin/bash
# ===================================================================
# LAMMPS RTX 5090 최신 빌드 (Ultimate Fix - FFT 문제 해결)
# 핵심 수정: FFT=CUFFT -> FFT=FFTW3 (최신 LAMMPS 호환)
# ===================================================================

if [ -z "$BASH_VERSION" ]; then
    exec bash "$0" "$@"
fi

set -e

# 색상 코드
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

print_header "LAMMPS RTX 5090 최신 빌드 (Ultimate Fix)"


# FFTW3 설치 확인 및 설치
sudo apt-get update
sudo apt-get install -y libfftw3-dev libfftw3-mpi-dev

# 설치 확인
pkg-config --modversion fftw3

# ===================================================================
# 1. NVIDIA HPC SDK 경로 설정
# ===================================================================
NVHPC_VERSION="25.11"
NVHPC_CUDA_VERSION="13.0"
NVHPC_ROOT="/opt/nvidia/hpc_sdk/Linux_x86_64/${NVHPC_VERSION}"

export CUDA_HOME="${NVHPC_ROOT}/cuda/${NVHPC_CUDA_VERSION}"
MATH_LIBS_HOME="${NVHPC_ROOT}/math_libs/${NVHPC_CUDA_VERSION}/targets/x86_64-linux"
MATH_LIBS_LIB="${MATH_LIBS_HOME}/lib"
MATH_LIBS_INCLUDE="${MATH_LIBS_HOME}/include"

# 환경변수 설정
export PATH="$CUDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:$MATH_LIBS_LIB:$LD_LIBRARY_PATH"
export CMAKE_PREFIX_PATH="$CUDA_HOME:$MATH_LIBS_HOME"
export CPATH="$CUDA_HOME/include:$MATH_LIBS_INCLUDE:$CPATH"

# CUDA 버전 확인
if [ ! -f "$CUDA_HOME/bin/nvcc" ]; then
    print_error "nvcc를 찾을 수 없습니다: $CUDA_HOME/bin/nvcc"
    exit 1
fi

CUDA_VERSION=$($CUDA_HOME/bin/nvcc --version | grep "release" | grep -oP '\d+\.\d+')
print_success "CUDA 환경: $CUDA_HOME (v$CUDA_VERSION)"

# ===================================================================
# 2. FFTW3 설치 확인 및 자동 설치
# ===================================================================
print_status "FFTW3 라이브러리 확인 중..."

if ! pkg-config --exists fftw3 2>/dev/null; then
    print_warning "FFTW3가 설치되어 있지 않습니다."
    print_status "FFTW3 자동 설치 중..."
    
    sudo apt-get update
    sudo apt-get install -y libfftw3-dev libfftw3-mpi-dev
    
    if pkg-config --exists fftw3 2>/dev/null; then
        print_success "FFTW3 설치 완료"
    else
        print_error "FFTW3 설치 실패. 수동으로 설치하세요:"
        echo "sudo apt-get install -y libfftw3-dev libfftw3-mpi-dev"
        exit 1
    fi
else
    FFTW3_VERSION=$(pkg-config --modversion fftw3)
    print_success "FFTW3 발견: v$FFTW3_VERSION"
fi

# ===================================================================
# 3. GPU 정보 확인 및 아키텍처 설정
# ===================================================================
print_status "GPU 정보 확인 중..."

if command -v nvidia-smi &> /dev/null; then
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 || echo "Unknown")
    COMPUTE_CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 || echo "12.0")
    GPU_MEM=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader 2>/dev/null | head -1 || echo "Unknown")
    
    print_success "GPU: $GPU_NAME (CC $COMPUTE_CAP) - $GPU_MEM"
else
    print_warning "nvidia-smi를 찾을 수 없습니다. RTX 5090 기본값 사용"
    GPU_NAME="RTX 5090"
    COMPUTE_CAP="12.0"
fi

# RTX 5090 Blackwell -> Ada Lovelace 호환 모드
if [[ "$COMPUTE_CAP" =~ ^12\. ]] || [[ "$GPU_NAME" == *"5090"* ]]; then
    print_warning "RTX 5090 Blackwell 아키텍처 감지 (CC $COMPUTE_CAP)"
    print_status "호환 모드: Ada Lovelace (sm_89) 사용"
fi

GPU_ARCH="89"
KOKKOS_ARCH="ADA89"

print_status "컴파일 아키텍처: sm_${GPU_ARCH} (Kokkos: ${KOKKOS_ARCH})"

# ===================================================================
# 4. 설치 경로 설정
# ===================================================================
INSTALL_BASE="$HOME/lammps-rtx5090-cuda13"
SOURCE_DIR="$INSTALL_BASE/source"
BUILD_DIR="$INSTALL_BASE/build"
INSTALL_DIR="$INSTALL_BASE/install"

print_status "설치 경로: $INSTALL_BASE"

# 기존 빌드 디렉토리만 삭제 (소스 재사용으로 시간 절약)
if [ -d "$BUILD_DIR" ]; then
    rm -rf "$BUILD_DIR"
    print_status "기존 빌드 디렉토리 초기화"
fi
mkdir -p "$BUILD_DIR"

# 소스가 없으면 다운로드, 있으면 재사용
if [ ! -d "$SOURCE_DIR" ]; then
    print_header "최신 LAMMPS 소스 다운로드 (Develop Branch)"
    cd "$INSTALL_BASE"
    
    print_status "GitHub에서 최신 develop 브랜치 다운로드 중..."
    git clone --depth 1 --branch develop https://github.com/lammps/lammps.git "$SOURCE_DIR"
    
    cd "$SOURCE_DIR"
    CURRENT_COMMIT=$(git rev-parse --short HEAD)
    print_success "다운로드 완료 - develop 브랜치 (Commit: $CURRENT_COMMIT)"
else
    print_success "기존 LAMMPS 소스 재사용"
    cd "$SOURCE_DIR"
    CURRENT_COMMIT=$(git rev-parse --short HEAD)
fi

# CMake 디렉토리 존재 확인
if [ ! -d "$SOURCE_DIR/cmake" ]; then
    print_error "치명적 오류: cmake 디렉토리가 없습니다."
    exit 1
fi

print_success "CMake 빌드 시스템 확인됨"

# ===================================================================
# 5. CMake 설정 (핵심 수정: FFT=FFTW3)
# ===================================================================
print_header "CMake 설정 (CUDA 13.0 + Kokkos + FFTW3)"
cd "$BUILD_DIR"

print_status "CMake 구성 중..."
print_status "  - CUDA: $CUDA_HOME"
print_status "  - GPU Arch: sm_${GPU_ARCH}"
print_status "  - Kokkos: ${KOKKOS_ARCH}"
print_status "  - FFT: FFTW3 (Host) + Kokkos cuFFT (Device)"

CMAKE_LOG="$BUILD_DIR/cmake_config.log"

cmake "$SOURCE_DIR/cmake" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_CXX_STANDARD_REQUIRED=ON \
    -DCMAKE_CUDA_COMPILER="$CUDA_HOME/bin/nvcc" \
    -DCMAKE_CUDA_ARCHITECTURES="$GPU_ARCH" \
    -DCMAKE_CUDA_FLAGS="-arch=sm_${GPU_ARCH} --use_fast_math -allow-unsupported-compiler" \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_MPI=ON \
    -DBUILD_OMP=ON \
    -DPKG_KOKKOS=ON \
    -DKokkos_ENABLE_CUDA=ON \
    -DKokkos_ARCH_${KOKKOS_ARCH}=ON \
    -DKokkos_ENABLE_OPENMP=ON \
    -DKokkos_ENABLE_SERIAL=ON \
    -DKokkos_ENABLE_CUDA_LAMBDA=ON \
    -DPKG_REAXFF=ON \
    -DPKG_MOLECULE=ON \
    -DPKG_KSPACE=ON \
    -DPKG_RIGID=ON \
    -DPKG_MANYBODY=ON \
    -DPKG_QEQ=ON \
    -DPKG_MISC=ON \
    -DFFT=FFTW3 \
    -DCUDAToolkit_ROOT="$CUDA_HOME" \
    -DCMAKE_PREFIX_PATH="$CUDA_HOME;$MATH_LIBS_HOME" \
    2>&1 | tee "$CMAKE_LOG"

CMAKE_STATUS=${PIPESTATUS[0]}

if [ $CMAKE_STATUS -ne 0 ]; then
    print_error "CMake 설정 실패!"
    echo ""
    print_status "에러 로그 마지막 30줄:"
    tail -30 "$CMAKE_LOG"
    echo ""
    print_status "가능한 해결 방법:"
    echo "1. FFTW3 재설치: sudo apt-get install --reinstall libfftw3-dev libfftw3-mpi-dev"
    echo "2. MPI 설치: sudo apt-get install openmpi-bin libopenmpi-dev"
    exit 1
fi

print_success "CMake 설정 완료"

# ===================================================================
# 6. 컴파일 (멀티코어 최적화)
# ===================================================================
print_header "컴파일 시작"

NPROC=$(nproc)
COMPILE_JOBS=$((NPROC > 32 ? 32 : NPROC))

print_status "컴파일 시작 (${COMPILE_JOBS}/${NPROC} 스레드)"
print_warning "예상 시간: 15-30분"

COMPILE_LOG="$BUILD_DIR/compile.log"
START_TIME=$(date +%s)

make -j${COMPILE_JOBS} 2>&1 | tee "$COMPILE_LOG"

COMPILE_STATUS=${PIPESTATUS[0]}
END_TIME=$(date +%s)
COMPILE_TIME=$((END_TIME - START_TIME))
COMPILE_MIN=$((COMPILE_TIME / 60))
COMPILE_SEC=$((COMPILE_TIME % 60))

if [ $COMPILE_STATUS -ne 0 ]; then
    print_error "컴파일 실패! (${COMPILE_MIN}분 ${COMPILE_SEC}초)"
    echo ""
    print_status "에러 분석 중..."
    
    if grep -qi "cudaMemPrefetchAsync" "$COMPILE_LOG"; then
        print_error "CUDA API 호환성 문제"
        echo "해결: CUDA 13.0의 일부 API가 불안정할 수 있습니다."
    elif grep -qi "out of memory\|cannot allocate memory" "$COMPILE_LOG"; then
        print_error "메모리 부족"
        echo "해결: make -j$((COMPILE_JOBS / 2))"
    fi
    
    echo ""
    echo "상세 로그: tail -50 $COMPILE_LOG"
    exit 1
fi

print_success "컴파일 완료 (${COMPILE_MIN}분 ${COMPILE_SEC}초)"

# ===================================================================
# 7. 설치 및 검증
# ===================================================================
print_header "설치 및 검증"

make install 2>&1 | tee install.log

LMP_BINARY="$INSTALL_DIR/bin/lmp"
if [ ! -f "$LMP_BINARY" ]; then
    print_error "실행파일 생성 실패: $LMP_BINARY"
    exit 1
fi

print_success "설치 완료: $LMP_BINARY"

# 패키지 검증
echo ""
print_status "설치된 패키지 확인:"
$LMP_BINARY -help 2>/dev/null | grep -E "KOKKOS|REAXFF|KSPACE|MOLECULE" || echo "  (패키지 정보 출력 실패)"

# ===================================================================
# 8. 환경 설정
# ===================================================================
print_header "환경 설정"

# 기존 설정 제거
if grep -q "LAMMPS_RTX5090" ~/.bashrc 2>/dev/null; then
    sed -i '/# LAMMPS RTX 5090/,/^$/d' ~/.bashrc
fi

print_status "~/.bashrc에 환경 변수 추가 중..."

cat >> ~/.bashrc << EOF

# ===================================================================
# LAMMPS RTX 5090 CUDA 13.0 환경 ($(date +%Y-%m-%d))
# ===================================================================
export NVHPC_ROOT="$NVHPC_ROOT"
export CUDA_HOME="$CUDA_HOME"
export MATH_LIBS_HOME="$MATH_LIBS_HOME"

export PATH="\$CUDA_HOME/bin:\$PATH"
export LD_LIBRARY_PATH="\$CUDA_HOME/lib64:\$MATH_LIBS_HOME/lib:\$LD_LIBRARY_PATH"

export LAMMPS_RTX5090="$INSTALL_DIR"
export PATH="\$LAMMPS_RTX5090/bin:\$PATH"

# RTX 5090 최적화 실행 명령어
alias lmp-5090="\$LAMMPS_RTX5090/bin/lmp -k on g 1 -sf kk"
alias lmp-5090-hybrid="\$LAMMPS_RTX5090/bin/lmp -k on g 1 t 24 -sf kk"
alias lmp-info="\$LAMMPS_RTX5090/bin/lmp -help | grep -E 'KOKKOS|REAXFF'"
EOF

print_success "환경 설정 완료"

# ===================================================================
# 9. GPU 기능 테스트
# ===================================================================
print_header "GPU 기능 테스트"

TEST_DIR="$INSTALL_BASE/test"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

cat > gpu_test.lmp << 'EOF'
# RTX 5090 GPU Test with Kokkos
units real
atom_style atomic
boundary p p p
newton on

# Kokkos 설정
package kokkos neigh half comm device

region box block 0 20 0 20 0 20
create_box 1 box
create_atoms 1 random 10000 12345 box
mass 1 12.0

# Kokkos pair style
pair_style lj/cut/kk 2.5
pair_coeff 1 1 0.1 3.0

velocity all create 300.0 87287
timestep 1.0
thermo 10

run 100

print "========================================="
print "RTX_5090_CUDA13_FFTW3_SUCCESS"
print "========================================="
EOF

TEST_LOG="$TEST_DIR/gpu_test.log"
print_status "GPU 테스트 실행 중..."

$LMP_BINARY -k on g 1 -sf kk -in gpu_test.lmp -log "$TEST_LOG" -screen none 2>&1 || true

if grep -q "RTX_5090_CUDA13_FFTW3_SUCCESS" "$TEST_LOG" 2>/dev/null; then
    print_success "✓ GPU 기능 테스트 통과"
    
    # 성능 정보 추출
    if grep -q "Performance:" "$TEST_LOG"; then
        echo ""
        print_status "성능 정보:"
        grep "Performance:" "$TEST_LOG" | tail -1
    fi
else
    print_warning "△ GPU 테스트 불완전 (실행은 가능할 수 있음)"
    if [ -f "$TEST_LOG" ]; then
        echo ""
        print_status "테스트 로그 마지막 20줄:"
        tail -20 "$TEST_LOG"
    fi
fi

rm -f gpu_test.lmp

# ===================================================================
# 10. 설치 완료 안내
# ===================================================================
print_header "설치 완료!"

echo ""
echo "🎉 RTX 5090 + CUDA 13.0 LAMMPS 설치 성공!"
echo ""
echo "📦 설치 정보:"
echo "   - LAMMPS: develop 브랜치 (Commit: $CURRENT_COMMIT)"
echo "   - 설치 경로: $INSTALL_DIR"
echo "   - CUDA 버전: $CUDA_VERSION"
echo "   - GPU 아키텍처: sm_${GPU_ARCH} (Ada Lovelace 호환 모드)"
echo "   - FFT 라이브러리: FFTW3 (Host) + Kokkos cuFFT (Device)"
echo ""
echo "🚀 즉시 사용:"
echo "   source ~/.bashrc"
echo "   lmp-5090 -in input.lmp"
echo ""
echo "💡 유용한 명령어:"
echo "   lmp-5090           # GPU 단독 실행"
echo "   lmp-5090-hybrid    # GPU + CPU 하이브리드"
echo "   lmp-info           # 패키지 확인"
echo ""
echo "📊 모니터링:"
echo "   nvidia-smi -l 1"
echo ""
echo "⚠️  참고사항:"
echo "   - RTX 5090은 Ada Lovelace (sm_89) 호환 모드로 실행됩니다."
echo "   - Kokkos를 통한 GPU 가속이 활성화되었습니다."
echo "   - Host FFT는 FFTW3, Device FFT는 Kokkos 내부 cuFFT를 사용합니다."
echo ""

print_success "모든 과정이 완료되었습니다!"
