import re
from tqdm import tqdm

def uniquify_pdb_atom_names(input_pdb, output_pdb):
    """
    PDB 원자명을 중복 없이 유니크하게 변환 (C001, O002 등)
    ATOM/HETATM 라인에만 적용
    """
    atom_counts = {}

    with open(input_pdb, 'r') as fin, open(output_pdb, 'w') as fout:
        for line in tqdm(fin, desc="Processing PDB"):
            if line.startswith(("ATOM", "HETATM")):
                # Extract element (77-78 column) or guess from atom name
                element = line[76:78].strip()
                if not element:
                    element = re.sub('[^A-Za-z]', '', line[12:16]).strip()
                    element = element[0] if element else "X"

                atom_counts[element] = atom_counts.get(element, 0) + 1
                new_atom_name = f"{element}{atom_counts[element]:03d}"  # C001 형식
                new_atom_name = new_atom_name[:4].ljust(4)  # PDB column width 유지

                # Replace atom name columns (13-16)
                new_line = line[:12] + new_atom_name + line[16:]
                fout.write(new_line)
            else:
                fout.write(line)

    print(f"[OK] Unique atom PDB saved → {output_pdb}")


if __name__ == "__main__":
    uniquify_pdb_atom_names("ionpairs.pdb", "ionpairs_unique.pdb")

