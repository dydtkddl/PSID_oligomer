#!/bin/bash
# run_bcmbp_research.sh - 순차 실행 스크립트

LAMMPS="/home/yongsang/downloads/lammps/build-reaxff/lmp"

echo "=========================================="
echo "BCMBP 연구용 시뮬레이션 실행"
echo "=========================================="

# 1단계: 벌크 평형화 (필수)
echo "1단계: 벌크 평형화 시작..."
mpirun -np 10 $LAMMPS -in 01_bulk_equilibration.lmp -log 01_equilibration.log

if [ $? -eq 0 ]; then
    echo "✅ 벌크 평형화 완료"
else
    echo "❌ 벌크 평형화 실패"
    exit 1
fi

# 2단계: 연구 시나리오 선택
echo ""
echo "연구 시나리오 선택:"
echo "1) 열분해 반응"
echo "2) 기계적 물성"
echo "3) 유리전이온도"
echo "4) 모든 시나리오"
read -p "선택 (1-4): " choice

case $choice in
    1)
        echo "열분해 시뮬레이션 시작..."
        mpirun -np 10 $LAMMPS -in 02_thermal_decomposition.lmp -log 02_thermal.log
        ;;
    2)
        echo "기계적 물성 측정 시작..."
        mpirun -np 10 $LAMMPS -in 03_mechanical_properties.lmp -log 03_mechanical.log
        ;;
    3)
        echo "Tg 계산 시작..."
        mpirun -np 10 $LAMMPS -in 04_glass_transition.lmp -log 04_tg.log
        ;;
    4)
        echo "모든 시나리오 실행..."
        mpirun -np 10 $LAMMPS -in 02_thermal_decomposition.lmp -log 02_thermal.log &
        mpirun -np 10 $LAMMPS -in 03_mechanical_properties.lmp -log 03_mechanical.log &
        mpirun -np 10 $LAMMPS -in 04_glass_transition.lmp -log 04_tg.log &
        wait
        ;;
esac

echo "=========================================="
echo "✅ 시뮬레이션 완료!"
echo "결과 분석을 시작하세요."
echo "=========================================="
