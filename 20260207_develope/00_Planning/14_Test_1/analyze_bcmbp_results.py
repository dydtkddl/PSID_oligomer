#!/usr/bin/env python3
# analyze_bcmbp_results.py - 결과 자동 분석

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def analyze_stress_strain():
    """응력-변형률 곡선 분석"""
    if Path("stress_strain.dat").exists():
        print("📊 응력-변형률 분석")
        data = np.loadtxt("stress_strain.dat")
        strain = data[:, 0]
        stress = data[:, 1]
        
        # 영률 계산 (초기 기울기)
        linear_region = strain < 0.02  # 2% 이하 선형 구간
        if np.sum(linear_region) > 10:
            young_modulus = np.polyfit(strain[linear_region], 
                                     stress[linear_region], 1)[0]
            print(f"  영률: {young_modulus:.2f} GPa")
        
        # 그래프 저장
        plt.figure(figsize=(10, 6))
        plt.plot(strain * 100, stress, 'b-', linewidth=2)
        plt.xlabel('변형률 (%)')
        plt.ylabel('응력 (GPa)')
        plt.title('BCMBP 응력-변형률 곡선')
        plt.grid(True, alpha=0.3)
        plt.savefig('stress_strain_curve.png', dpi=300, bbox_inches='tight')
        print("  그래프 저장: stress_strain_curve.png")

def analyze_glass_transition():
    """유리전이온도 분석"""
    if Path("density_temperature.dat").exists():
        print("📊 유리전이온도 분석")
        data = np.loadtxt("density_temperature.dat")
        temp = data[:, 0]
        density = data[:, 1]
        
        # 그래프 저장
        plt.figure(figsize=(10, 6))
        plt.plot(temp, density, 'ro-', linewidth=2, markersize=6)
        plt.xlabel('온도 (K)')
        plt.ylabel('밀도 (g/cm³)')
        plt.title('BCMBP 밀도-온도 관계 (Tg 계산용)')
        plt.grid(True, alpha=0.3)
        plt.savefig('density_temperature.png', dpi=300, bbox_inches='tight')
        print("  그래프 저장: density_temperature.png")
        print("  Tg 계산: 기울기 변화점을 찾아 수동 결정")

if __name__ == "__main__":
    print("=== BCMBP 시뮬레이션 결과 분석 ===\n")
    
    analyze_stress_strain()
    print()
    analyze_glass_transition()
    
    print("\n=== 분석 완료 ===")
