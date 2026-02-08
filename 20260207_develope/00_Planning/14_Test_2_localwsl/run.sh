#!/bin/bash
# ===================================================================
# BCMBP 연구용 시뮬레이션 통합 실행 스크립트 (01~05)
# GPU 자동 감지, 에러 처리, 결과 정리 포함
# ===================================================================

set -e  # 에러 발생 시 스크립트 중단

# 색상 정의
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

print_header() { echo -e "${CYAN}========================================${NC}"; echo -e "${CYAN}$1${NC}"; echo -e "${CYAN}========================================${NC}"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_info() { echo -e "${BLUE}[INFO]${NC} $1"; }

# ===================================================================
# 1. 실행 환경 자동 감지
# ===================================================================
print_header "BCMBP 연구용 시뮬레이션 실행"

# GPU LAMMPS 경로 확인 (우선순위 순)
GPU_PATHS=(
    "$HOME/lammps-gpu-reaxff/install/bin/lmp"
    "$HOME/lammps_cufft_fixed/install/bin/lmp"
    "$(which lmp-gpu 2>/dev/null || echo '')"
)

CPU_LAMMPS="/home/yongsang/downloads/lammps/build-reaxff/lmp"

# GPU 실행 환경 설정
USE_GPU=0
LAMMPS_CMD=""

for gpu_path in "${GPU_PATHS[@]}"; do
    if [ -n "$gpu_path" ] && [ -f "$gpu_path" ] && command -v nvidia-smi &> /dev/null; then
        LAMMPS_CMD="$gpu_path -k on g 1 -sf kk"
        USE_GPU=1
        print_success "GPU 가속 모드: $gpu_path"
        break
    fi
done

# GPU를 찾지 못한 경우 CPU 사용
if [ $USE_GPU -eq 0 ]; then
    if [ -f "$CPU_LAMMPS" ]; then
        LAMMPS_CMD="mpirun -np 10 $CPU_LAMMPS"
        print_warning "CPU 모드 사용: $CPU_LAMMPS"
    else
        print_error "LAMMPS 실행 파일을 찾을 수 없습니다!"
        exit 1
    fi
fi

print_info "실행 명령어: $LAMMPS_CMD"

# GPU 상태 확인 (GPU 모드인 경우)
if [ $USE_GPU -eq 1 ]; then
    echo ""
    print_info "GPU 상태:"
    nvidia-smi --query-gpu=name,temperature.gpu,memory.used,memory.total --format=csv,noheader,nounits | head -1
fi

# ===================================================================
# 2. 필수 파일 확인
# ===================================================================
echo ""
print_header "필수 파일 확인"

REQUIRED_FILES=("bcmbp.data" "11_CHOCl.lammps.ff")
missing_files=0

for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "$file" ]; then
        print_success "✓ $file"
    else
        print_error "✗ $file 없음!"
        missing_files=$((missing_files + 1))
    fi
done

if [ $missing_files -gt 0 ]; then
    print_error "$missing_files개의 필수 파일이 없습니다."
    exit 1
fi

# 05_oligomerization.lmp 자동 생성 (없는 경우)
if [ ! -f "05_oligomerization.lmp" ]; then
    print_warning "05_oligomerization.lmp 없음. 자동 생성 중..."
    
    cat > 05_oligomerization.lmp << 'EOF'
# ===================================================================
# BCMBP 올리고머 형성 시뮬레이션
# 전략: 단계적 가열 → 반응 유도 → 냉각 안정화
# ===================================================================

units           real
atom_style      full
boundary        p p p

print "=========================================="
print "BCMBP 올리고머 형성 시뮬레이션"
print "메커니즘: C-Cl 결합 파괴 → 라디칼 결합"
print "=========================================="

# 1. 평형 벌크 구조 로드
read_data       bcmbp_bulk_equilibrated.data

# 2. ReaxFF 설정
pair_style      reaxff NULL safezone 4.0 mincap 300
pair_coeff      * * 11_CHOCl.lammps.ff C H Cl
fix             qeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff

neighbor        2.5 bin
neigh_modify    every 5 delay 0 check yes

timestep        0.2  # 반응 시뮬레이션용 작은 타임스텝

# 3. 화학 반응 분석 설정
fix             species all reaxff/species 100 1000 1000 species_oligo.out element C H Cl
fix             bonds all reaxff/bonds 1000 bonds_oligo.reax

# 4. 출력 설정
thermo          1000
thermo_style    custom step temp press pe ke etotal vol density

# ===================================================================
# 단계 1: C-Cl 결합 활성화 (300K → 600K, 10 ps)
# ===================================================================
print ""
print "=== 단계 1: C-Cl 결합 활성화 (300K → 600K) ==="
velocity        all create 300.0 123456 dist gaussian

fix             heat1 all nvt temp 300.0 600.0 100.0
dump            d1 all custom 2000 05_activation.lammpstrj id type xu yu zu q
dump_modify     d1 sort id

run             50000  # 10 ps
undump          d1
unfix           heat1

# ===================================================================
# 단계 2: 라디칼 형성 (600K → 1000K, 30 ps)
# ===================================================================
print ""
print "=== 단계 2: 라디칼 형성 (600K → 1000K) ==="

# 추가 가열
fix             heat2 all nvt temp 600.0 1000.0 100.0
run             25000  # 5 ps
unfix           heat2

# 고온 반응 구간
print "고온 반응 구간 (1000K, 25 ps)"
fix             react all nvt temp 1000.0 1000.0 100.0
dump            d2 all custom 1000 05_reaction.lammpstrj id type xu yu zu q
dump_modify     d2 sort id

run             125000  # 25 ps
undump          d2
unfix           react

# ===================================================================
# 단계 3: 올리고머 형성 (1000K → 500K, 25 ps)
# ===================================================================
print ""
print "=== 단계 3: 올리고머 형성 (냉각) ==="

fix             cool1 all nvt temp 1000.0 500.0 100.0
dump            d3 all custom 1000 05_oligomerization.lammpstrj id type xu yu zu q
dump_modify     d3 sort id

run             125000  # 25 ps
undump          d3
unfix           cool1

# ===================================================================
# 단계 4: 안정화 (500K → 300K, 15 ps)
# ===================================================================
print ""
print "=== 단계 4: 안정화 (500K → 300K) ==="

fix             cool2 all nvt temp 500.0 300.0 100.0
run             50000  # 10 ps
unfix           cool2

# 최종 평형화
fix             final all nvt temp 300.0 300.0 100.0
run             25000  # 5 ps
unfix           final

# ===================================================================
# 결과 저장
# ===================================================================
unfix           species
unfix           bonds

write_data      bcmbp_oligomer_final.data

print ""
print "=========================================="
print "✅ 올리고머 형성 시뮬레이션 완료!"
print "총 시뮬레이션 시간: 100 ps"
print ""
print "핵심 분석 파일:"
print "  - species_oligo.out: 화학종 변화"
print "  - bonds_oligo.reax: 결합 생성/파괴"
print "  - 05_oligomerization.lammpstrj: 형성 과정"
print "=========================================="
EOF

    print_success "05_oligomerization.lmp 생성 완료"
fi

# ===================================================================
# 3. 벌크 평형화 확인 (1단계)
# ===================================================================
echo ""
print_header "1단계: 벌크 평형화"

if [ ! -f "bcmbp_bulk_equilibrated.data" ]; then
    print_info "벌크 평형화 시작..."
    
    start_time=$(date +%s)
    $LAMMPS_CMD -in 01_bulk_equilibration.lmp -log 01_equilibration.log
    end_time=$(date +%s)
    runtime=$((end_time - start_time))
    
    if [ $? -eq 0 ]; then
        print_success "벌크 평형화 완료 (${runtime}초)"
    else
        print_error "벌크 평형화 실패!"
        echo "로그 확인: tail -20 01_equilibration.log"
        exit 1
    fi
else
    print_success "기존 평형 구조 발견. 벌크 평형화 스킵."
fi

# ===================================================================
# 4. 연구 시나리오 선택
# ===================================================================
echo ""
print_header "연구 시나리오 선택"

echo ""
echo "  1) 🔥 열분해 반응 (Thermal Decomposition)"
echo "  2) 🔨 기계적 물성 (Mechanical Properties)"  
echo "  3) ❄️  유리전이온도 (Glass Transition)"
echo "  4) 🔗 올리고머 형성 (Oligomerization) ⭐ 추천"
echo "  5) 🚀 모든 시나리오 순차 실행"
echo "  6) ⚡ 모든 시나리오 병렬 실행 (고성능 시스템)"
echo ""
read -p "선택 (1-6): " choice

# 시뮬레이션 실행 함수
run_simulation() {
    local input_file=$1
    local log_file=$2
    local description=$3
    
    echo ""
    print_info "$description 시작..."
    
    local start_time=$(date +%s)
    $LAMMPS_CMD -in "$input_file" -log "$log_file"
    local exit_code=$?
    local end_time=$(date +%s)
    local runtime=$((end_time - start_time))
    
    if [ $exit_code -eq 0 ]; then
        print_success "$description 완료 (${runtime}초)"
        return 0
    else
        print_error "$description 실패!"
        echo "로그 확인: tail -20 $log_file"
        return 1
    fi
}

# ===================================================================
# 5. 선택된 시나리오 실행
# ===================================================================
case $choice in
    1)
        run_simulation "02_thermal_decomposition.lmp" "02_thermal.log" "열분해 반응"
        ;;
    2)
        run_simulation "03_mechanical_properties.lmp" "03_mechanical.log" "기계적 물성 측정"
        ;;
    3)
        run_simulation "04_glass_transition.lmp" "04_tg.log" "유리전이온도 계산"
        ;;
    4)
        run_simulation "05_oligomerization.lmp" "05_oligo.log" "올리고머 형성"
        ;;
    5)
        print_header "모든 시나리오 순차 실행"
        
        run_simulation "02_thermal_decomposition.lmp" "02_thermal.log" "2/4: 열분해 반응"
        run_simulation "03_mechanical_properties.lmp" "03_mechanical.log" "3/4: 기계적 물성"
        run_simulation "04_glass_transition.lmp" "04_tg.log" "4/4: 유리전이온도"
        run_simulation "05_oligomerization.lmp" "05_oligo.log" "5/5: 올리고머 형성"
        ;;
    6)
        print_header "모든 시나리오 병렬 실행"
        print_warning "시스템 리소스를 많이 사용합니다!"
        
        $LAMMPS_CMD -in 02_thermal_decomposition.lmp -log 02_thermal.log &
        PID1=$!
        $LAMMPS_CMD -in 03_mechanical_properties.lmp -log 03_mechanical.log &
        PID2=$!
        $LAMMPS_CMD -in 04_glass_transition.lmp -log 04_tg.log &
        PID3=$!
        $LAMMPS_CMD -in 05_oligomerization.lmp -log 05_oligo.log &
        PID4=$!
        
        print_info "백그라운드 실행 중... (PID: $PID1, $PID2, $PID3, $PID4)"
        wait $PID1 $PID2 $PID3 $PID4
        print_success "모든 병렬 시뮬레이션 완료"
        ;;
    *)
        print_error "잘못된 선택: $choice"
        exit 1
        ;;
esac

# ===================================================================
# 6. 결과 정리 및 다음 단계 안내
# ===================================================================
echo ""
print_header "시뮬레이션 완료"

# 결과 파일 확인
echo ""
print_info "생성된 주요 파일:"
ls -lh *.lammpstrj *.out *.reax *.dat 2>/dev/null | head -10 | while read line; do
    echo "  $line"
done

echo ""
print_header "다음 단계: 결과 분석"
echo ""
echo "📊 **분석 방법:**"
echo ""
echo "1. 올리고머 형성 확인:"
echo "   tail -20 species_oligo.out"
echo ""
echo "2. 궤적 시각화 (OVITO):"
echo "   ovito 05_oligomerization.lammpstrj"
echo ""
echo "3. GPU 사용률 확인:"
if [ $USE_GPU -eq 1 ]; then
    echo "   nvidia-smi"
else
    echo "   (CPU 모드로 실행됨)"
fi
echo ""
echo "4. 결과 분석 스크립트:"
echo "   python3 analyze_bcmbp_results.py"
echo ""

print_success "모든 작업이 완료되었습니다! 🎉"

