# ReaxFF 파라미터 실제 사용 및 MD 소프트웨어 호환성 테스트 계획

검증된 `11_CHOCl.ff` 파라미터를 실제 MD 소프트웨어에서 돌려보는 **"Hello World"** 수준의 테스트 계획을 제시합니다. 목표는 **30분 내에 완료 가능한 간단한 테스트**로 소프트웨어 호환성을 확인하는 것입니다.

## **1단계: 소프트웨어 호환성 평가**

### **ReaxFF 지원 MD 소프트웨어 현황**

| 소프트웨어  | ReaxFF 지원 | 설치 난이도 | 속도      | 추천도     | 비고                     |
| ----------- | ----------- | ----------- | --------- | ---------- | ------------------------ |
| **LAMMPS**  | ✅ 완벽     | 중간        | 빠름      | ⭐⭐⭐⭐⭐ | 오픈소스, 가장 널리 사용 |
| **AMS**     | ✅ 완벽     | 쉬움        | 매우 빠름 | ⭐⭐⭐⭐   | 상용, trial 가능         |
| **GULP**    | ✅ 제한적   | 어려움      | 중간      | ⭐⭐⭐     | 오픈소스, 고체 특화      |
| **GROMACS** | ❌ 미지원   | -           | -         | ❌         | ReaxFF 없음              |
| **NAMD**    | ❌ 미지원   | -           | -         | ❌         | ReaxFF 없음              |

**결론: LAMMPS를 1차 테스트 대상으로, AMS를 2차 검증 대상으로 선정**

## **2단계: 테스트 시스템 설계**

### **벤질클로라이드 (C₇H₇Cl) 선택 근거**

**BCMBP와의 직접적 연관성:**

- BCMBP 구조: `ClCH₂-Ar-Ar-CH₂Cl`
- 벤질클로라이드: `ClCH₂-Ar` (BCMBP의 핵심 단위)
- C-Cl 결합 파라미터 직접 테스트 가능
- 작은 시스템 (15원자)으로 빠른 검증

**테스트 시나리오:**

```
1. 에너지 최소화 (구조 안정성 확인)
2. 300K NVT 평형화 (5 ps) - 기본 동역학
3. 고온 테스트 1000K (2 ps) - C-Cl 결합 활성화 확인
```

## **3단계: LAMMPS 빠른 테스트 (20분 완료)**

### **3-1. 벤질클로라이드 구조 파일 생성**

다음 내용을 `benzyl_chloride.data`로 저장:

```lammps
LAMMPS data file - Benzyl Chloride (C7H7Cl) ReaxFF Test

15 atoms
4 atom types

-10.0 10.0 xlo xhi
-10.0 10.0 ylo yhi
-10.0 10.0 zlo zhi

Masses

1 12.0000  # C
2 1.0080   # H
3 15.9990  # O (사용 안함)
4 35.4500  # Cl

Atoms  # full

# 벤젠 고리 탄소 (6개)
1  1 1  0.0  -1.185  0.374  0.000  # C1 (벤젠)
2  1 1  0.0  -0.439  1.602  0.000  # C2
3  1 1  0.0   0.957  1.666  0.000  # C3
4  1 1  0.0   1.608  0.437  0.000  # C4
5  1 1  0.0   0.862 -0.791  0.000  # C5
6  1 1  0.0  -0.534 -0.823  0.000  # C6

# 메틸렌 탄소
7  1 1  0.0  -2.700  0.280  0.000  # CH2

# 수소 원자 (7개)
8  1 2  0.0  -0.908  2.528  0.000  # H (벤젠)
9  1 2  0.0   1.503  2.604  0.000  # H
10 1 2  0.0   2.693  0.487  0.000  # H
11 1 2  0.0   1.362 -1.728  0.000  # H
12 1 2  0.0  -1.077 -1.761  0.000  # H
13 1 2  0.0  -3.090  1.299  0.000  # H (메틸렌)
14 1 2  0.0  -3.100 -0.700  0.000  # H (메틸렌)

# 염소 원자
15 1 4  0.0  -3.500 -0.800  0.000  # Cl
```

### **3-2. LAMMPS 입력 스크립트**

다음 내용을 `test_reaxff.lmp`로 저장:

```lammps
# ===================================================================
# 11_CHOCl.ff ReaxFF 파라미터 호환성 테스트
# 시스템: 벤질클로라이드 (C7H7Cl)
# ===================================================================

# 기본 설정
units real
atom_style charge
boundary p p p

# 구조 읽기
read_data benzyl_chloride.data

# ===================================================================
# ReaxFF 설정 - 핵심 부분!
# ===================================================================
pair_style reaxff NULL safezone 3.0 mincap 150
pair_coeff * * 11_CHOCl.ff C H O Cl

# 전하 평형 (ReaxFF 필수)
fix qeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff

# 출력 설정
thermo 50
thermo_style custom step temp pe ke etotal press vol

# ===================================================================
# 테스트 1: 에너지 최소화
# ===================================================================
print "=== 에너지 최소화 시작 ==="
minimize 1.0e-4 1.0e-6 1000 10000
print "=== 에너지 최소화 완료 ==="

# ===================================================================
# 테스트 2: 300K 평형화 (5 ps)
# ===================================================================
reset_timestep 0
timestep 0.25  # fs

velocity all create 300.0 12345 dist gaussian
fix nvt1 all nvt temp 300.0 300.0 100.0

dump 1 all custom 100 benzyl_300K.lammpstrj id type x y z q
dump_modify 1 sort id

print "=== 300K 평형화 시작 ==="
run 20000  # 5 ps
print "=== 300K 평형화 완료 ==="

unfix nvt1
undump 1

# ===================================================================
# 테스트 3: 고온 안정성 (1000K, 2 ps)
# ===================================================================
fix nvt2 all nvt temp 300.0 1000.0 100.0  # 가열
run 2000  # 0.5 ps 가열

unfix nvt2
fix nvt3 all nvt temp 1000.0 1000.0 100.0

dump 2 all custom 100 benzyl_1000K.lammpstrj id type x y z q
dump_modify 2 sort id

print "=== 1000K 고온 테스트 시작 ==="
run 8000  # 2 ps at 1000K
print "=== 1000K 고온 테스트 완료 ==="

unfix nvt3
undump 2

# 최종 구조 저장
write_data benzyl_final.data

print "======================================"
print "✅ ReaxFF 호환성 테스트 완료!"
print "파라미터: 11_CHOCl.ff"
print "시스템: 벤질클로라이드 (15 원자)"
print "총 시뮬레이션: 7.5 ps"
print "======================================"
```

### **3-3. 실행 및 결과 확인**

```bash
# 실행 (5-10분 소요)
lammps -in test_reaxff.lmp -log reaxff_test.log

# 또는 병렬 실행
mpirun -np 4 lammps -in test_reaxff.lmp -log reaxff_test.log
```

**성공 판별 기준:**

```bash
# 1. 에러 메시지 확인
grep -i "error" reaxff_test.log
# 출력 없으면 성공

# 2. 완료 메시지 확인
grep "완료" reaxff_test.log
# 3개 단계 모두 완료 메시지 있어야 함

# 3. 에너지 안정성 확인
tail -20 reaxff_test.log
# 마지막 온도가 1000K 근처, 에너지가 유한값

# 4. 출력 파일 생성 확인
ls -lh benzyl_*.lammpstrj benzyl_final.data
```

**예상 정상 출력:**

```
Step Temp PotEng KinEng TotEng Press
...
20000 299.5 -1250.3 49.8 -1200.5 12.3  # 300K 단계
...
28000 998.7 -1180.2 166.1 -1014.1 45.6  # 1000K 단계
✅ ReaxFF 호환성 테스트 완료!
```

## **4단계: AMS 검증 테스트 (선택사항, 10분)**

LAMMPS 테스트가 성공하면 AMS에서도 확인:

```python
#!/usr/bin/env python3
"""AMS ReaxFF 호환성 빠른 테스트"""

from scm.plams import *

def quick_ams_test():
    init()

    # 벤질클로라이드 분자 (SMILES)
    mol = from_smiles('ClCC1=CC=CC=C1')

    # ReaxFF 설정
    s = Settings()
    s.input.ReaxFF.ForceField = '11_CHOCl.ff'
    s.input.ams.Task = 'MolecularDynamics'
    s.input.ams.MolecularDynamics.NSteps = 2000  # 0.5 ps
    s.input.ams.MolecularDynamics.TimeStep = 0.25
    s.input.ams.MolecularDynamics.InitialVelocities.Temperature = 300
    s.input.ams.MolecularDynamics.Thermostat.Type = 'Berendsen'
    s.input.ams.MolecularDynamics.Thermostat.Temperature = 300

    # 실행
    job = AMSJob(molecule=mol, settings=s, name='ams_reaxff_test')
    result = job.run()

    if result.ok():
        print("✅ AMS ReaxFF 호환성 확인!")
        print(f"평균 온도: {result.get_history_property('Temperature')[-1]:.1f} K")
        return True
    else:
        print("❌ AMS 테스트 실패")
        return False

    finish()

if __name__ == "__main__":
    quick_ams_test()
```

## **5단계: 문제 해결 가이드**

### **일반적인 오류와 해결책**

| 오류 메시지                         | 원인               | 해결 방법                                   |
| ----------------------------------- | ------------------ | ------------------------------------------- |
| `Cannot open ReaxFF potential file` | 파일 경로 문제     | 절대 경로 사용: `/full/path/to/11_CHOCl.ff` |
| `Illegal pair_style command`        | ReaxFF 패키지 없음 | LAMMPS 재컴파일: `make yes-reaxff`          |
| `atoms lost` (시뮬레이션 폭발)      | timestep 너무 큼   | `timestep 0.1`로 감소                       |
| `NaN` 에너지                        | 파라미터 문제      | 파일 형식 재확인                            |

### **성능 벤치마크**

**예상 실행 시간:**

- 에너지 최소화: 1-2분
- 300K MD (5 ps): 2-3분
- 1000K MD (2 ps): 1-2분
- **총 소요 시간: 5-10분**

## **6단계: 다음 단계 로드맵**

### **호환성 테스트 성공 시**

```
✅ 11_CHOCl.ff → LAMMPS 호환성 확인
✅ C-Cl 결합 파라미터 정상 작동
✅ 고온에서 구조 안정성 확인
➡️ BCMBP 단량체 테스트 준비
➡️ 소규모 가교중합 시뮬레이션 설계
```

### **실패 시 대안 전략**

```
❌ 파라미터 파일 문제 → 형식 재검토
❌ LAMMPS 설치 문제 → Docker 이미지 사용
❌ 구조 불안정 → 더 간단한 분자(CH₃Cl)로 시작
```

## **최종 실행 체크리스트**

```bash
# 준비 사항 (2분)
[  ] LAMMPS 설치 확인: lammps --version
[  ] 11_CHOCl.ff 파일 위치 확인: ls -l 11_CHOCl.ff
[  ] 작업 디렉토리 생성: mkdir reaxff_test && cd reaxff_test

# 파일 생성 (3분)
[  ] benzyl_chloride.data 생성
[  ] test_reaxff.lmp 생성

# 실행 및 확인 (10분)
[  ] LAMMPS 실행: lammps -in test_reaxff.lmp
[  ] 로그 파일 에러 확인: grep -i error *.log
[  ] 출력 파일 확인: ls *.lammpstrj *.data
[  ] 온도/에너지 안정성 확인

# 총 소요 시간: 15-20분
```

이 테스트가 성공하면 `11_CHOCl.ff` 파라미터는 **실제 MD 시뮬레이션에서 완벽하게 작동**하며, BCMBP 가교중합 연구로 안전하게 진행할 수 있습니다!
