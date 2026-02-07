#!/usr/bin/env python3
"""
RDKit으로 정확한 BCMBP 구조 생성
"""
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def create_bcmbp_data():
    # BCMBP SMILES
    smiles = "ClCc1ccc(cc1)c2ccc(cc2)CCl"
    
    # 분자 생성 및 최적화
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    # 3D 좌표 생성 (ETKDG 방법 - 더 현실적인 구조)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
    
    conf = mol.GetConformer()
    n_atoms = mol.GetNumAtoms()
    
    print(f"BCMBP 분자 생성 완료")
    print(f"분자식: C14H12Cl2")
    print(f"총 원자 수: {n_atoms}")
    
    # 원자 타입 매핑
    type_map = {"C": 1, "H": 2, "O": 3, "Cl": 4}
    
    atoms_data = []
    type_counts = {"C": 0, "H": 0, "Cl": 0}
    
    for i in range(n_atoms):
        atom = mol.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        symbol = atom.GetSymbol()
        
        atoms_data.append({
            'id': i + 1,
            'symbol': symbol,
            'type': type_map[symbol],
            'coords': [pos.x, pos.y, pos.z]
        })
        
        type_counts[symbol] += 1
    
    print(f"\n원자 구성:")
    print(f"  C (탄소): {type_counts['C']}개")
    print(f"  H (수소): {type_counts['H']}개")
    print(f"  Cl (염소): {type_counts['Cl']}개")
    
    # 박스 크기 계산
    coords = np.array([atom['coords'] for atom in atoms_data])
    margin = 10.0
    box_min = coords.min(axis=0) - margin
    box_max = coords.max(axis=0) + margin
    
    # LAMMPS data 파일 작성
    with open("bcmbp.data", "w") as f:
        f.write("LAMMPS data file - BCMBP 4,4-bis(chloromethyl)-1,1-biphenyl\n\n")
        f.write(f"{n_atoms} atoms\n")
        f.write("4 atom types\n\n")
        f.write(f"{box_min[0]:.3f} {box_max[0]:.3f} xlo xhi\n")
        f.write(f"{box_min[1]:.3f} {box_max[1]:.3f} ylo yhi\n")
        f.write(f"{box_min[2]:.3f} {box_max[2]:.3f} zlo zhi\n\n")
        
        f.write("Masses\n\n")
        f.write("1 12.0000\n")   # C
        f.write("2 1.0080\n")    # H
        f.write("3 15.9990\n")   # O (ReaxFF 호환성)
        f.write("4 35.4500\n\n") # Cl
        
        f.write("Atoms\n\n")
        for atom in atoms_data:
            mol_id = 1
            charge = 0.0  # ReaxFF가 계산
            x, y, z = atom['coords']
            f.write(f"{atom['id']:3d} {mol_id} {atom['type']} {charge:.4f} "
                   f"{x:8.4f} {y:8.4f} {z:8.4f}\n")
    
    print(f"\n✅ bcmbp.data 파일 생성 완료!")
    return True

if __name__ == "__main__":
    create_bcmbp_data()
