import os
import argparse
import logging
import numpy as np
from tqdm import tqdm
from ase.io import read, write
from xtb.ase.calculator import XTB

def setup_logging(logfile):
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

def main(args):
    # ===== Load XYZ frames =====
    logging.info(f"Loading input XYZ: {args.input}")
    db = read(args.input, ":")
    logging.info(f"Loaded {len(db)} frames")

    # ===== Initialize xTB calculator =====
    xtb_calc = XTB(method=args.method)
    logging.info(f"xTB method initialized: {args.method}")

    energies, forces = [], []

    # ===== xTB calculations =====
    for i, at in enumerate(tqdm(db, desc="Running xTB single-point")):
        try:
            at.calc = xtb_calc
            e = at.get_potential_energy()
            f = at.get_forces()

            at.info["energy_xtb"] = e
            at.arrays["forces_xtb"] = f

            energies.append(e)
            forces.append(f)
            logging.info(f"Frame {i}: energy = {e:.6f} eV, forces OK")
        except Exception as err:
            logging.error(f"Frame {i} failed: {err}")

    energies = np.array(energies)
    forces = np.array(forces)

    # ===== Save extended XYZ =====
    write(args.out_xyz, db)
    logging.info(f"Saved extended XYZ: {args.out_xyz}")

    # ===== Save NPZ dataset =====
    np.savez(args.out_npz, energies=energies, forces=forces)
    logging.info(f"Saved dataset NPZ: {args.out_npz}")

    print("\n✅ xTB dataset generation complete!")
    print(f"📦 Extended XYZ → {args.out_xyz}")
    print(f"💾 NumPy dataset → {args.out_npz}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build xTB dataset from multi-frame XYZ")
    parser.add_argument("--input", type=str, default="xtbopt/all_ionpairs.xyz",
                        help="Input multi-frame XYZ file")
    parser.add_argument("--out_xyz", type=str, default="xtbopt/xtb_dataset.xyz",
                        help="Output extended XYZ file")
    parser.add_argument("--out_npz", type=str, default="xtbopt/xtb_dataset.npz",
                        help="Output NPZ dataset")
    parser.add_argument("--method", type=str, default="GFN2-xTB",
                        help="xTB method (GFN2-xTB, GFN1-xTB, GFN-FF)")
    parser.add_argument("--log", type=str, default="xtb_dataset.log",
                        help="Log file name")

    args = parser.parse_args()
    setup_logging(args.log)
    main(args)

