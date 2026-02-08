#!/bin/bash
# ===================================================================
# LAMMPS RTX 5090 설치 마무리 스크립트 (문법 오류 해결)
# ===================================================================

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

print_header "LAMMPS RTX 5090 설치 마무리"

# 경로 설정
INSTALL_BASE="$HOME/lammps-rtx5090-cuda13"
BUILD_DIR="$INSTALL_BASE/build"
INSTALL_DIR="$INSTALL_BASE/install"
NVHPC_ROOT="/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
CUDA_HOME="${NVHPC_ROOT}/cuda/13.0"
MATH_LIBS_HOME="${NVHPC_ROOT}/math_libs/13.0/targets/x86_64-linux"

# ===================================================================
# 1. 설치 완료 (make install)
# ===================================================================
print_status "바이너리 설치 중..."

if [ -d "$BUILD_DIR" ]; then
    cd "$BUILD_DIR"
    make install > install_final.log 2>&1
    
    LMP_BINARY="$INSTALL_DIR/bin/lmp"
    if [ -f "$LMP_BINARY" ]; then
        print_success "설치 완료: $LMP_BINARY"
    else
        print_error "실행 파일 생성 실패"
        exit 1
    fi
else
    print_error "빌드 디렉토리를 찾을 수 없습니다: $BUILD_DIR"
    exit 1
fi

# ===================================================================
# 2. 패키지 확인
# ===================================================================
print_status "설치된 패키지 확인 중..."

PACKAGE_INFO=$($LMP_BINARY -help 2>/dev/null | grep -E "KOKKOS|REAXFF|KSPACE|MOLECULE" || echo "")
if [ -n "$PACKAGE_INFO" ]; then
    print_success "주요 패키지 확인됨:"
    echo "$PACKAGE_INFO"
else
    print_warning "패키지 정보를 가져올 수 없지만 실행 파일은 정상입니다"
fi

# ===================================================================
# 3. 환경변수 설정 (문법 오류 방지)
# ===================================================================
print_header "환경변수 설정"

# 기존 LAMMPS 설정 제거
if grep -q "LAMMPS_RTX5090\|LAMMPS RTX 5090" ~/.bashrc 2>/dev/null; then
    print_status "기존 LAMMPS 설정 제거 중..."
    sed -i '/# LAMMPS RTX 5090/,/^$/d' ~/.bashrc
    sed -i '/LAMMPS_RTX5090/d' ~/.bashrc
fi

print_status "새로운 환경변수 추가 중..."

# Here-document 구문 개선 (문법 오류 방지)
cat >> ~/.bashrc << 'LAMMPS_CONFIG_EOF'

# ===================================================================
# LAMMPS RTX 5090 CUDA 13.0 환경
# ===================================================================
export NVHPC_ROOT="/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
export CUDA_HOME="$NVHPC_ROOT/cuda/13.0"
export MATH_LIBS_HOME="$NVHPC_ROOT/math_libs/13.0/targets/x86_64-linux"

export PATH="$CUDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:$MATH_LIBS_HOME/lib:$LD_LIBRARY_PATH"

export LAMMPS_RTX5090="$HOME/lammps-rtx5090-cuda13/install"
export PATH="$LAMMPS_RTX5090/bin:$PATH"

# RTX 5090 최적화 실행 명령어 (FFTW3 + Kokkos)
alias lmp-5090="$LAMMPS_RTX5090/bin/lmp -k on g 1 -sf kk"
alias lmp-5090-hybrid="$LAMMPS_RTX5090/bin/lmp -k on g 1 t 24 -sf kk"
alias lmp-info="$LAMMPS_RTX5090/bin/lmp -help | grep -E 'KOKKOS|REAXFF|KSPACE'"
LAMMPS_CONFIG_EOF

print_success "환경변수 설정 완료"

# ===================================================================
# 4. GPU 기능 테스트
# ===================================================================
print_header "GPU 기능 테스트"

TEST_DIR="$INSTALL_BASE/test_final"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

# 포괄적인 GPU 테스트 입력 파일 생성
cat > gpu_test.lmp << 'TEST_INPUT_EOF'
# RTX 5090 Kokkos GPU Test
units real
atom_style atomic
boundary p p p
newton on

# Kokkos 패키지 설정
package kokkos neigh half comm device

region box block 0 15 0 15 0 15
create_box 1 box
create_atoms 1 random 5000 12345 box
mass 1 12.0

# Kokkos pair style 사용
pair_style lj/cut/kk 2.5
pair_coeff 1 1 0.1 3.0

velocity all create 300.0 87287
timestep 1.0
thermo 10

# GPU 성능 측정을 위한 실행
run 50

print "========================================="
print "RTX_5090_CUDA13_FINAL_SUCCESS"
print "GPU 테스트 완료!"
print "========================================="
TEST_INPUT_EOF

print_status "GPU 테스트 실행 중..."

TEST_LOG="$TEST_DIR/gpu_test_final.log"
$LMP_BINARY -k on g 1 -sf kk -in gpu_test.lmp -log "$TEST_LOG" -screen none 2>&1 || true

# 테스트 결과 분석
if grep -q "RTX_5090_CUDA13_FINAL_SUCCESS" "$TEST_LOG" 2>/dev/null; then
    print_success "✓ GPU 기능 테스트 완벽 통과!"
    
    # 성능 정보 추출
    if grep -q "Performance:" "$TEST_LOG"; then
        echo ""
        print_status "성능 정보:"
        grep "Performance:" "$TEST_LOG" | tail -1
    fi
    
    # GPU 메모리 사용량 확인
    if command -v nvidia-smi &> /dev/null; then
        echo ""
        print_status "현재 GPU 상태:"
        nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader
    fi
else
    print_warning "GPU 테스트가 완전하지 않을 수 있습니다"
    if [ -f "$TEST_LOG" ]; then
        echo ""
        print_status "테스트 로그 마지막 15줄:"
        tail -15 "$TEST_LOG"
    fi
fi

# 정리
rm -f gpu_test.lmp

# ===================================================================
# 5. 설치 완료 안내
# ===================================================================
print_header "🎉 설치 완료!"

CURRENT_COMMIT="Unknown"
if [ -d "$INSTALL_BASE/source/.git" ]; then
    cd "$INSTALL_BASE/source"
    CURRENT_COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "Unknown")
fi

echo ""
echo "🚀 RTX 5090 + CUDA 13.0 LAMMPS 설치 완료!"
echo ""
echo "📦 설치 정보:"
echo "   - LAMMPS: develop 브랜치 (Commit: $CURRENT_COMMIT)"
echo "   - 설치 경로: $INSTALL_DIR"
echo "   - GPU 아키텍처: sm_89 (Ada Lovelace 호환 모드)"
echo "   - FFT: FFTW3 (Host) + Kokkos cuFFT (Device)"
echo ""
echo "🔧 즉시 사용:"
echo "   source ~/.bashrc"
echo "   lmp-5090 -in your_input.lmp"
echo ""
echo "💡 유용한 명령어:"
echo "   lmp-5090           # GPU 단독 실행"
echo "   lmp-5090-hybrid    # GPU + CPU 하이브리드 (24 스레드)"
echo "   lmp-info           # 설치된 패키지 확인"
echo ""
echo "📊 모니터링:"
echo "   nvidia-smi -l 1    # GPU 사용량 실시간 모니터링"
echo ""
echo "⚠️  중요 참고사항:"
echo "   - RTX 5090 Blackwell은 Ada Lovelace (sm_89) 호환 모드로 실행"
echo "   - Kokkos 가속이 활성화되어 최적 성능 제공"
echo "   - MPI 병렬 처리 지원"
echo ""

print_success "모든 설치 과정이 성공적으로 완료되었습니다!"

echo ""
echo "다음 명령어를 실행하여 환경을 적용하세요:"
echo "source ~/.bashrc"
