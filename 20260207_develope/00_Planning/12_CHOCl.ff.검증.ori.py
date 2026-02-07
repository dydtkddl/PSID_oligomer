
#!/usr/bin/env python3
"""
ReaxFF 파라미터 파일 종합 검증 스크립트
11_CHOCl.ff 파일의 형식, 구조, BCMBP 적합성을 검증합니다.
"""

import os
import sys

def parse_clean_line(line):
    """주석 제거 및 토큰 추출"""
    clean = line.split('!')[0].strip()
    return clean.split() if clean else []

def validate_reaxff_file(filename):
    """
    ReaxFF 파라미터 파일 종합 검증

    Returns:
    --------
    tuple: (success: bool, info: dict)
    """

    print(f"{'='*70}")
    print(f"ReaxFF 파라미터 파일 검증: {filename}")
    print(f"{'='*70}")

    # 파일 존재 확인
    if not os.path.exists(filename):
        print(f"❌ 파일을 찾을 수 없습니다: {filename}")
        return False, {}

    file_info = {
        'filename': filename,
        'file_size': os.path.getsize(filename),
        'n_general': 0,
        'n_atoms': 0,
        'atoms': [],
        'n_bonds': 0,
        'n_offdiag': 0,
        'n_angles': 0,
        'n_torsions': 0,
        'n_hbonds': 0,
        'has_chlorine': False,
        'has_ccl_bond': False
    }

    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        print(f"✓ 파일 읽기 성공 ({len(lines)} 라인, {file_info['file_size']} bytes)")

        # 데이터 라인만 추출
        data_lines = []
        for i, line in enumerate(lines):
            tokens = parse_clean_line(line)
            if tokens:
                data_lines.append((i+1, tokens))  # (라인번호, 토큰들)

        print(f"✓ 데이터 라인 추출 완료 ({len(data_lines)}개)\n")

        current_idx = 0

        # 1. 헤더 확인
        if lines:
            header = lines[0].strip()
            print(f"[헤더] {header}")
            if not data_lines[0][1][0].isdigit():
                current_idx = 1  # 첫 라인이 헤더인 경우

        # 2. 일반 파라미터
        print(f"\n[섹션 1] 일반 파라미터")
        line_num, tokens = data_lines[current_idx]
        file_info['n_general'] = int(float(tokens[0]))
        print(f"  개수: {file_info['n_general']} (라인 {line_num})")
        current_idx += 1

        # 일반 파라미터 값들 읽기
        params_count = 0
        while params_count < file_info['n_general']:
            if current_idx >= len(data_lines):
                raise ValueError(f"일반 파라미터 {params_count}/{file_info['n_general']}에서 파일 종료")

            line_num, tokens = data_lines[current_idx]
            try:
                # 각 라인의 첫 번째 값만 파라미터로 계산 (ReaxFF 표준)
                float(tokens[0])
                params_count += 1
                current_idx += 1
            except (ValueError, IndexError):
                raise ValueError(f"라인 {line_num}에서 일반 파라미터 파싱 실패: {tokens}")

        print(f"  ✓ {params_count}개 파라미터 검증 완료")

        # 3. 원자 타입
        print(f"\n[섹션 2] 원자 타입")
        line_num, tokens = data_lines[current_idx]
        file_info['n_atoms'] = int(float(tokens[0]))
        print(f"  개수: {file_info['n_atoms']} (라인 {line_num})")
        current_idx += 1

        # 각 원자는 4줄로 구성 (심볼 + 3줄의 파라미터)
        for atom_idx in range(file_info['n_atoms']):
            if current_idx >= len(data_lines):
                raise ValueError(f"원자 {atom_idx+1}에서 파일 종료")

            # 원자 심볼 라인
            line_num, tokens = data_lines[current_idx]
            atom_symbol = tokens[0]
            file_info['atoms'].append(atom_symbol)
            print(f"  - 원자 {atom_idx+1}: {atom_symbol} (라인 {line_num})")

            # 원자당 4줄 (심볼 라인 포함) 소비
            current_idx += 4

            if atom_symbol == 'Cl':
                file_info['has_chlorine'] = True

        if file_info['has_chlorine']:
            print(f"  ✅ Cl(염소) 원자 발견 - BCMBP 시뮬레이션 호환")
        else:
            print(f"  ❌ Cl 원자 없음 - BCMBP에 부적합")

        # 4. 결합 파라미터
        print(f"\n[섹션 3] 결합 파라미터")
        line_num, tokens = data_lines[current_idx]
        file_info['n_bonds'] = int(float(tokens[0]))
        print(f"  개수: {file_info['n_bonds']} (라인 {line_num})")
        current_idx += 1

        # 각 결합은 2줄로 구성
        ccl_bonds = []
        for bond_idx in range(file_info['n_bonds']):
            if current_idx >= len(data_lines):
                raise ValueError(f"결합 {bond_idx+1}에서 파일 종료")

            # 결합 정의 라인 (at1 at2 De_sigma ...)
            line_num, tokens = data_lines[current_idx]
            try:
                at1, at2 = int(tokens[0]), int(tokens[1])

                # 원자 심볼 매핑
                atom1_symbol = file_info['atoms'][at1-1] if at1 <= len(file_info['atoms']) else '?'
                atom2_symbol = file_info['atoms'][at2-1] if at2 <= len(file_info['atoms']) else '?'
                bond_name = f"{atom1_symbol}-{atom2_symbol}"

                print(f"  - 결합 {bond_idx+1}: {bond_name} ({at1}-{at2})")

                # C-Cl 결합 확인
                if (atom1_symbol == 'C' and atom2_symbol == 'Cl') or \
                   (atom1_symbol == 'Cl' and atom2_symbol == 'C'):
                    file_info['has_ccl_bond'] = True
                    ccl_bonds.append(bond_name)

            except (ValueError, IndexError):
                print(f"  ⚠ 결합 {bond_idx+1} 파싱 실패 (라인 {line_num})")

            current_idx += 2  # 결합당 2줄 소비

        if file_info['has_ccl_bond']:
            print(f"  ✅ C-Cl 결합 파라미터 발견: {ccl_bonds}")
        else:
            print(f"  ❌ C-Cl 결합 파라미터 없음")

        # 5-7. 나머지 섹션들 (간단 검증)
        sections = [
            ('Off-diagonal 항', 'n_offdiag', 1),
            ('각도 파라미터', 'n_angles', 1),
            ('비틀림 파라미터', 'n_torsions', 1),
            ('수소결합 파라미터', 'n_hbonds', 1)
        ]

        for section_num, (name, key, lines_per_item) in enumerate(sections, 4):
            print(f"\n[섹션 {section_num}] {name}")

            if current_idx >= len(data_lines):
                print(f"  ⚠ 섹션 없음 (파일 종료)")
                file_info[key] = 0
                continue

            line_num, tokens = data_lines[current_idx]
            count = int(float(tokens[0]))
            file_info[key] = count
            print(f"  개수: {count} (라인 {line_num})")
            current_idx += 1

            # 해당 개수만큼 라인 건너뛰기
            current_idx += count * lines_per_item
            print(f"  ✓ {count}개 항목 건너뛰기 완료")

        # 최종 검증
        print(f"\n{'='*70}")
        print(f"검증 결과 요약")
        print(f"{'='*70}")

        # 예상값과 비교
        expected = {
            'n_general': 41,
            'n_atoms': 4,
            'atoms': ['C', 'H', 'O', 'Cl'],
            'n_bonds': 10,
            'n_offdiag': 6,
            'n_angles': 29,
            'n_torsions': 32,
            'n_hbonds': 2
        }

        all_good = True

        for key, exp_val in expected.items():
            if key == 'atoms':
                actual = file_info[key]
                match = actual == exp_val
                status = "✅" if match else "⚠️"
                print(f"{status} 원자 타입: {actual} {'(정상)' if match else f'(예상: {exp_val})'}")
            else:
                actual = file_info[key]
                match = actual == exp_val
                status = "✅" if match else "⚠️"
                print(f"{status} {key}: {actual} {'(정상)' if match else f'(예상: {exp_val})'}")

            if not match:
                all_good = False

        # BCMBP 적합성 확인
        print(f"\n[BCMBP 적합성 검사]")
        bcmbp_ready = file_info['has_chlorine'] and file_info['has_ccl_bond']
        status = "✅" if bcmbp_ready else "❌"
        print(f"{status} Cl 원자: {'있음' if file_info['has_chlorine'] else '없음'}")
        print(f"{status} C-Cl 결합: {'있음' if file_info['has_ccl_bond'] else '없음'}")

        final_success = all_good and bcmbp_ready

        print(f"\n{'='*70}")
        if final_success:
            print("🎉 검증 성공! BCMBP 시뮬레이션에 사용 가능합니다.")
        else:
            print("⚠️ 일부 문제 발견. 파일을 재확인하세요.")
        print(f"{'='*70}")

        return final_success, file_info

    except Exception as e:
        print(f"\n❌ 검증 중 오류 발생:")
        print(f"   {str(e)}")
        import traceback
        traceback.print_exc()
        return False, file_info


def quick_fix_suggestions(filename, info):
    """파일 수정 제안"""
    print(f"\n[수정 제안]")

    if info.get('n_atoms') != 4:
        print("- 원자 개수가 4개가 아닙니다. C, H, O, Cl 확인 필요")

    if not info.get('has_chlorine'):
        print("- Cl 원자가 없습니다. 원자 섹션에 'Cl' 추가 필요")

    if not info.get('has_ccl_bond'):
        print("- C-Cl 결합 파라미터가 없습니다. 결합 섹션 확인 필요")

    print("- PDF 복사 시 줄바꿈 오류가 있을 수 있습니다.")
    print("- 숫자가 합쳐지거나 분리되었는지 확인하세요.")


if __name__ == "__main__":
    # 파일명 설정
    filename = "11_CHOCl.ff"

    # 명령행 인자 처리
    if len(sys.argv) > 1:
        filename = sys.argv[1]

    # 검증 실행
    success, info = validate_reaxff_file(filename)

    # 실패 시 수정 제안
    if not success:
        quick_fix_suggestions(filename, info)

    # 종료 코드 설정
    sys.exit(0 if success else 1)
