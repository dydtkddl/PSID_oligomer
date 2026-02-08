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
