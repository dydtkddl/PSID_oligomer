#!/usr/bin/env python3
"""
BCMBP ReaxFF 시뮬레이션 완전 분석 및 시각화
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import warnings
warnings.filterwarnings('ignore')

# 한글 폰트 설정 (선택사항)
try:
    plt.rcParams['font.family'] = ['DejaVu Sans', 'Arial']
except:
    pass

print("="*60)
print("BCMBP ReaxFF 시뮬레이션 종합 분석")
print("="*60)

# ===================================================================
# 1. 에너지 분석
# ===================================================================
print("\n[1/5] 에너지 분석 중...")

try:
    data = np.loadtxt('energy_vs_time.dat', skiprows=1)
    steps = data[:, 0]
    temps = data[:, 1]
    pe = data[:, 2]
    ke = data[:, 3]
    etot = data[:, 4]
    press = data[:, 5]

    time_ps = steps * 0.25 / 1000  # fs to ps

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # 온도 프로파일
    axes[0, 0].plot(time_ps, temps, 'b-', linewidth=1.5)
    axes[0, 0].axhline(300, color='r', linestyle='--', alpha=0.7, label='300K')
    axes[0, 0].axhline(800, color='orange', linestyle='--', alpha=0.7, label='800K')
    axes[0, 0].set_xlabel('Time (ps)')
    axes[0, 0].set_ylabel('Temperature (K)')
    axes[0, 0].set_title('Temperature Profile')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # 전체 에너지
    axes[0, 1].plot(time_ps, etot, 'g-', linewidth=1.5)
    axes[0, 1].set_xlabel('Time (ps)')
    axes[0, 1].set_ylabel('Total Energy (kcal/mol)')
    axes[0, 1].set_title('Total Energy Evolution')
    axes[0, 1].grid(True, alpha=0.3)

    # PE vs KE
    axes[1, 0].plot(time_ps, pe, 'r-', linewidth=1.5, label='Potential Energy')
    axes[1, 0].plot(time_ps, ke, 'b-', linewidth=1.5, label='Kinetic Energy')
    axes[1, 0].set_xlabel('Time (ps)')
    axes[1, 0].set_ylabel('Energy (kcal/mol)')
    axes[1, 0].set_title('Potential vs Kinetic Energy')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    # 압력
    axes[1, 1].plot(time_ps, press, 'purple', linewidth=1.5)
    axes[1, 1].axhline(0, color='k', linestyle='--', alpha=0.5)
    axes[1, 1].set_xlabel('Time (ps)')
    axes[1, 1].set_ylabel('Pressure (atm)')
    axes[1, 1].set_title('Pressure Evolution')
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('energy_analysis.png', dpi=300, bbox_inches='tight')
    print("✓ 저장: energy_analysis.png")

except Exception as e:
    print(f"⚠ 에너지 분석 건너뜀: {e}")

# ===================================================================
# 2. XYZ Trajectory 분석
# ===================================================================
print("\n[2/5] XYZ trajectory 분석 중...")

def read_xyz_trajectory(filename):
    """다중 프레임 XYZ 파일 읽기"""
    frames = []
    atoms_info = []
    
    try:
        with open(filename, 'r') as f:
            while True:
                line = f.readline().strip()
                if not line:
                    break
                    
                n_atoms = int(line)
                comment = f.readline().strip()
                
                coords = []
                elements = []
                for _ in range(n_atoms):
                    parts = f.readline().split()
                    elements.append(parts[0])
                    coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
                
                frames.append(np.array(coords))
                if not atoms_info:  # 첫 프레임에서만 원소 정보 저장
                    atoms_info = elements
                    
    except Exception as e:
        print(f"XYZ 파일 읽기 실패: {e}")
        return [], []
    
    return frames, atoms_info

try:
    frames, elements = read_xyz_trajectory('bcmbp_complete.xyz')
    
    if frames:
        print(f"✓ {len(frames)} 프레임 로드 완료")
        
        # 구조 변화 분석
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        
        # 첫 프레임
        coords_0 = frames[0]
        axes[0].scatter(coords_0[:, 0], coords_0[:, 1], c='blue', s=50, alpha=0.7)
        axes[0].set_xlabel('X (Å)')
        axes[0].set_ylabel('Y (Å)')
        axes[0].set_title('Initial Structure (t=0)')
        axes[0].set_aspect('equal')
        axes[0].grid(True, alpha=0.3)
        
        # 중간 프레임
        mid_idx = len(frames) // 2
        coords_mid = frames[mid_idx]
        axes[1].scatter(coords_mid[:, 0], coords_mid[:, 1], c='green', s=50, alpha=0.7)
        axes[1].set_xlabel('X (Å)')
        axes[1].set_ylabel('Y (Å)')
        axes[1].set_title(f'Middle Structure (frame {mid_idx})')
        axes[1].set_aspect('equal')
        axes[1].grid(True, alpha=0.3)
        
        # 최종 프레임
        coords_final = frames[-1]
        axes[2].scatter(coords_final[:, 0], coords_final[:, 1], c='red', s=50, alpha=0.7)
        axes[2].set_xlabel('X (Å)')
        axes[2].set_ylabel('Y (Å)')
        axes[2].set_title('Final Structure')
        axes[2].set_aspect('equal')
        axes[2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('structure_evolution.png', dpi=300, bbox_inches='tight')
        print("✓ 저장: structure_evolution.png")
        
except Exception as e:
    print(f"⚠ XYZ 분석 건너뜀: {e}")

# ===================================================================
# 3. 결합 길이 분석
# ===================================================================
print("\n[3/5] 결합 길이 분석 중...")

try:
    bond_data = np.loadtxt('bond_lengths_800K.dat', skiprows=1)
    
    plt.figure(figsize=(14, 6))
    
    # 결합 길이 분포
    plt.subplot(1, 2, 1)
    all_bonds = bond_data[:, 1:].flatten()
    plt.hist(all_bonds, bins=50, alpha=0.7, edgecolor='black', color='skyblue')
    plt.xlabel('Bond Length (Å)')
    plt.ylabel('Frequency')
    plt.title('Bond Length Distribution at 800K')
    plt.grid(True, alpha=0.3)
    
    # 평균 결합 길이 시간 변화
    plt.subplot(1, 2, 2)
    mean_bonds = np.mean(bond_data[:, 1:], axis=1)
    time_steps = bond_data[:, 0] * 0.25 / 1000  # ps 변환
    plt.plot(time_steps, mean_bonds, 'r-', linewidth=2)
    plt.xlabel('Time (ps)')
    plt.ylabel('Average Bond Length (Å)')
    plt.title('Average Bond Length Evolution (800K)')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('bond_analysis.png', dpi=300, bbox_inches='tight')
    print("✓ 저장: bond_analysis.png")
    
except Exception as e:
    print(f"⚠ 결합 분석 건너뜀: {e}")

# ===================================================================
# 4. 웹 기반 인터랙티브 시각화 (선택사항)
# ===================================================================
print("\n[4/5] 인터랙티브 시각화 준비 중...")

try:
    # py3Dmol 사용 가능 시 HTML 생성
    import py3Dmol
    
    def create_interactive_html():
        with open('bcmbp_complete.xyz', 'r') as f:
            xyz_content = f.read()
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>BCMBP ReaxFF Trajectory</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
</head>
<body>
    <div id="viewer" style="width: 800px; height: 600px;"></div>
    <div>
        <button onclick="viewer.animate({{loop: 'forward', interval: 100}})">Play</button>
        <button onclick="viewer.animate({{loop: 'stop'}})">Stop</button>
    </div>
    <script>
        let viewer = $3Dmol.createViewer('viewer');
        viewer.addModelsAsFrames(`{xyz_content}`, 'xyz');
        viewer.setStyle({{}}, {{stick: {{radius: 0.15}}, sphere: {{scale: 0.3}}}});
        viewer.zoomTo();
        viewer.render();
    </script>
</body>
</html>
        """
        
        with open('bcmbp_interactive.html', 'w') as f:
            f.write(html_content)
        print("✓ 저장: bcmbp_interactive.html")
    
    create_interactive_html()
    
except ImportError:
    print("⚠ py3Dmol 없음 - pip install py3Dmol로 설치 후 재실행")

# ===================================================================
# 5. 종합 리포트 생성
# ===================================================================
print("\n[5/5] 종합 리포트 생성 중...")

with open('analysis_summary.txt', 'w', encoding='utf-8') as f:
    f.write("="*60 + "\n")
    f.write("BCMBP ReaxFF 시뮬레이션 분석 종합 리포트\n")
    f.write("="*60 + "\n\n")
    
    try:
        f.write(f"총 시뮬레이션 시간: {time_ps[-1]:.2f} ps\n")
        f.write(f"총 스텝 수: {int(steps[-1])}\n")
        f.write(f"프레임 수: {len(frames)}\n\n")
        
        f.write("온도 통계:\n")
        f.write(f"  평균: {np.mean(temps):.2f} K\n")
        f.write(f"  최소: {np.min(temps):.2f} K\n")
        f.write(f"  최대: {np.max(temps):.2f} K\n\n")
        
        f.write("에너지 통계 (kcal/mol):\n")
        f.write(f"  PE 평균: {np.mean(pe):.2f}\n")
        f.write(f"  KE 평균: {np.mean(ke):.2f}\n")
        f.write(f"  Total 평균: {np.mean(etot):.2f}\n\n")
        
    except:
        f.write("데이터 분석 중 일부 오류 발생\n\n")
    
    f.write("생성된 파일 목록:\n")
    f.write("  📊 energy_analysis.png - 에너지 프로파일\n")
    f.write("  🔬 structure_evolution.png - 구조 변화\n")
    f.write("  🔗 bond_analysis.png - 결합 길이 분석\n")
    f.write("  🌐 bcmbp_interactive.html - 웹 시각화\n")
    f.write("  📄 analysis_summary.txt - 이 리포트\n")

print("✓ 저장: analysis_summary.txt")

print("\n" + "="*60)
print("✅ 완전 분석 완료!")
print("="*60)
print("🎯 시각화 방법:")
print("  VMD: vmd bcmbp_complete.dcd")
print("  OVITO: ovito bcmbp_complete.xyz")
print("  웹브라우저: bcmbp_interactive.html 열기")
print("  Python: 생성된 PNG 파일들 확인")
print("="*60)
