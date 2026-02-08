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
