#!/bin/bash
# BCMBP ReaxFF 테스트 자동 실행 스크립트

LAMMPS="/home/yongsang/downloads/lammps/build-reaxff/lmp"  # 또는 ~/downloads/lammps/build-cpu/lmp

echo "=========================================="
echo "BCMBP ReaxFF 테스트 자동 실행"
echo "=========================================="

# 1. LAMMPS 실행 가능 확인
if ! command -v $LAMMPS &> /dev/null; then
    echo "❌ LAMMPS를 찾을 수 없습니다"
    echo "PATH 확인 또는 절대경로 사용 필요"
    exit 1
fi
echo "✅ LAMMPS 실행파일 확인: $(which $LAMMPS)"

# 2. ReaxFF 패키지 확인
if $LAMMPS -help | grep -q "REAXFF"; then
    echo "✅ ReaxFF 패키지 설치 확인됨"
else
    echo "❌ ReaxFF 패키지가 설치되지 않았습니다"
    echo "LAMMPS 재컴파일이 필요합니다"
    exit 1
fi

# 3. 필수 파일 확인
echo ""
echo "=== 필수 파일 확인 ==="
missing_files=0
for file in test_bcmbp_reaxff.lmp bcmbp.data 11_CHOCl.lammps.ff; do
    if [ -f "$file" ]; then
        echo "✅ $file"
    else
        echo "❌ $file 없음!"
        missing_files=$((missing_files + 1))
    fi
done

if [ $missing_files -gt 0 ]; then
    echo "필수 파일이 없습니다. 파일을 생성한 후 다시 실행하세요."
    exit 1
fi

# 4. 테스트 실행
echo ""
echo "=== BCMBP 테스트 실행 ==="
echo "로그 파일: bcmbp_test.log"
echo "실시간 확인: tail -f bcmbp_test.log (별도 터미널에서)"
echo ""

$LAMMPS -in test_bcmbp_reaxff.lmp -log bcmbp_test.log

# 5. 결과 확인
echo ""
echo "=== 테스트 결과 ==="
if grep -q "완료" bcmbp_test.log 2>/dev/null; then
    echo "🎉 테스트 성공!"
    echo ""
    echo "생성된 파일:"
    ls -lh bcmbp_*.lammpstrj bcmbp_*.data 2>/dev/null || echo "출력 파일을 찾을 수 없습니다"
else
    echo "❌ 테스트 실패 또는 미완료"
    echo "로그 확인: tail -20 bcmbp_test.log"
fi

echo "=========================================="
