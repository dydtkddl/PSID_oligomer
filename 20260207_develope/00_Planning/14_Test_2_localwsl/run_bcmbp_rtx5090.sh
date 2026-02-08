#!/bin/bash
# ===================================================================
# RTX 5090 BCMBP: 실행 모드 선택 가능한 하이브리드 버전
# ===================================================================

echo "=========================================="
echo "BCMBP RTX 5090: 실행 모드 선택"
echo "01_평형화 + 05_올리고머화"
echo "=========================================="

# ===================================================================
# 실행 모드 선택
# ===================================================================
echo ""
echo "🔧 실행 모드 선택:"
echo "   1) GPU 전용 (RTX 5090만 사용)"
echo "   2) CPU+GPU 하이브리드 (권장: 최고 성능)"
echo "   3) 자동 선택 (시스템 최적화)"
echo ""
read -p "선택 [1-3, 기본값=2]: " MODE_CHOICE
MODE_CHOICE=${MODE_CHOICE:-2}

# CPU 정보 자동 감지
TOTAL_CORES=$(nproc)
PHYSICAL_CORES=$((TOTAL_CORES / 2))  # 하이퍼스레딩 고려

case $MODE_CHOICE in
    1)
        EXEC_MODE="GPU 전용"
        LAMMPS_ARGS="-k on g 1 -sf kk"
        echo "✅ GPU 전용 모드: RTX 5090만 사용"
        ;;
    2)
        EXEC_MODE="CPU+GPU 하이브리드"
        # 물리 코어의 75% 사용 (안정성 고려)
        CPU_THREADS=$((PHYSICAL_CORES * 3 / 4))
        [ $CPU_THREADS -lt 8 ] && CPU_THREADS=8
        [ $CPU_THREADS -gt 24 ] && CPU_THREADS=24
        
        LAMMPS_ARGS="-k on g 1 t $CPU_THREADS -sf kk"
        echo "✅ 하이브리드 모드: RTX 5090 + CPU ${CPU_THREADS}스레드"
        ;;
    3)
        EXEC_MODE="자동 최적화"
        # ReaxFF + RTX 5090 조합에서는 GPU 전용이 더 빠를 수 있음
        if [ $TOTAL_CORES -ge 32 ]; then
            CPU_THREADS=16
            LAMMPS_ARGS="-k on g 1 t $CPU_THREADS -sf kk"
            echo "✅ 자동 선택: 하이브리드 모드 (고성능 CPU 감지)"
        else
            LAMMPS_ARGS="-k on g 1 -sf kk"
            echo "✅ 자동 선택: GPU 전용 모드"
        fi
        ;;
esac

# ===================================================================
# 필수 파일 확인 (기존 코드와 동일)
# ===================================================================
REQUIRED_FILES=("bcmbp.data" "11_CHOCl.lammps.ff" "01_bulk_equilibration.lmp" "05_oligomerization.lmp")

echo ""
echo "📋 필수 파일 확인 중..."
for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "❌ 파일 없음: $file"
        exit 1
    else
        echo "✅ 파일 확인: $file"
    fi
done

# ===================================================================
# GPU 상태 확인
# ===================================================================
if command -v nvidia-smi &> /dev/null; then
    GPU_FREE=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null | head -1)
    echo ""
    echo "🖥️  GPU 여유 메모리: ${GPU_FREE} MB"
    echo "🔧 실행 모드: $EXEC_MODE"
fi

# ===================================================================
# 출력 디렉토리 및 실행
# ===================================================================
OUTPUT_DIR="results_$(echo $EXEC_MODE | tr ' ' '_')_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

LAMMPS_BINARY="$HOME/lammps-rtx5090-cuda13/install/bin/lmp"
LAMMPS_CMD="$LAMMPS_BINARY $LAMMPS_ARGS"

echo ""
echo "🚀 실행 명령어: $LAMMPS_CMD"
echo "📁 출력 디렉토리: $OUTPUT_DIR"
echo ""

# # ===================================================================
# # 단계별 실행 (기존 로직과 동일)
# # ===================================================================
# echo "=========================================="
# echo "🚀 단계 1/2: Bulk Equilibration 시작"
# echo "=========================================="

# START_TIME_01=$(date +%s)
# $LAMMPS_CMD -in 01_bulk_equilibration.lmp -log "$OUTPUT_DIR/01_equilibration.log"

# if [ $? -ne 0 ]; then
#     echo "❌ 단계 1 실패!"
#     exit 1
# fi

# END_TIME_01=$(date +%s)
# RUNTIME_01=$((END_TIME_01 - START_TIME_01))
# echo "✅ 단계 1 완료! (${RUNTIME_01}초)"

# # 중간 파일 복사
# cp *.data *.restart "$OUTPUT_DIR/" 2>/dev/null || true

echo ""
echo "=========================================="
echo "🚀 단계 2/2: Oligomerization 시작"
echo "=========================================="

START_TIME_05=$(date +%s)
$LAMMPS_CMD -in 05_oligomerization.lmp -log "$OUTPUT_DIR/05_oligomerization.log"

if [ $? -ne 0 ]; then
    echo "❌ 단계 2 실패!"
    exit 1
fi

END_TIME_05=$(date +%s)
RUNTIME_05=$((END_TIME_05 - START_TIME_05))
echo "✅ 단계 2 완료! (${RUNTIME_05}초)"

# ===================================================================
# 결과 정리
# ===================================================================
mv *.log *.lammpstrj *.dump *.reax *species*.out *.data *.restart "$OUTPUT_DIR/" 2>/dev/null || true

TOTAL_RUNTIME=$((RUNTIME_01 + RUNTIME_05))
TOTAL_MIN=$((TOTAL_RUNTIME / 60))
TOTAL_SEC=$((TOTAL_RUNTIME % 60))

echo ""
echo "=========================================="
echo "🎉 시뮬레이션 완료!"
echo "=========================================="
echo "⏱️  총 실행 시간: ${TOTAL_MIN}분 ${TOTAL_SEC}초"
echo "🔧 실행 모드: $EXEC_MODE"
echo "📁 결과 위치: $OUTPUT_DIR"
echo ""

# 성능 정보 추출 (가능한 경우)
if [ -f "$OUTPUT_DIR/05_oligomerization.log" ]; then
    echo "📊 성능 정보:"
    grep "Performance:" "$OUTPUT_DIR/05_oligomerization.log" | tail -1 || echo "   (성능 정보 없음)"
fi

echo ""
echo "✅ RTX 5090 BCMBP 시뮬레이션 완료! 🚀"
