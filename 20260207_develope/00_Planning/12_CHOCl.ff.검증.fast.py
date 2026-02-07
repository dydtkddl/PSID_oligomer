#!/usr/bin/env python3
"""
빠른 ReaxFF 파일 검증 (핵심만)
"""

def quick_check(filename="11_CHOCl.ff"):
    """5초 안에 끝나는 빠른 검증"""

    print(f"빠른 검증: {filename}")
    print("-" * 40)

    try:
        with open(filename, 'r') as f:
            lines = [line.split('!')[0].strip() for line in f.readlines() if line.split('!')[0].strip()]

        # 핵심 지표만 확인
        checks = []

        # 1. 일반 파라미터 (보통 2번째 라인)
        n_general = int(lines[0] if lines[0].isdigit() else lines[1])
        checks.append(("일반 파라미터", n_general, 41))

        # 2. Cl 원자 존재
        content = ' '.join(lines)
        has_cl = 'Cl ' in content and '35.45' in content
        checks.append(("Cl 원자", "있음" if has_cl else "없음", "있음"))

        # 3. C-Cl 결합 (1=C, 4=Cl)
        has_ccl = '1 4' in content
        checks.append(("C-Cl 결합", "있음" if has_ccl else "없음", "있음"))

        # 결과 출력
        all_ok = True
        for name, actual, expected in checks:
            ok = actual == expected
            status = "✅" if ok else "❌"
            print(f"{status} {name}: {actual}")
            if not ok:
                all_ok = False

        print("-" * 40)
        print(f"{'✅ 사용 가능' if all_ok else '❌ 문제 있음'}")

        return all_ok

    except Exception as e:
        print(f"❌ 오류: {e}")
        return False

if __name__ == "__main__":
    import sys
    filename = sys.argv[1] if len(sys.argv) > 1 else "11_CHOCl.ff"
    success = quick_check(filename)
    sys.exit(0 if success else 1)
