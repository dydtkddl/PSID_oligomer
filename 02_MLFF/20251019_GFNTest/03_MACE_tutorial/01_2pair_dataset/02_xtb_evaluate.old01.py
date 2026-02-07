import os
import argparse
import logging
from tqdm import tqdm
from ase.io import read, write
from xtb.ase.calculator import XTB
import numpy as np

def setup_logging(logfile):
    logging.basicConfig(
        filename=logfile,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        filemode='w'
    )

def run_xtb(args):
    # ASE로 XYZ 읽기
    logging.info(f"🔹 Loading XYZ frames from: {args.input}")
    frames = read(args.input, ":")

    # xTB 계산기 세팅
    logging.info(f"🔹 Initializing xTB method: {args.method}")
    xtb_calc = XTB(method=args.method)

    energies, forces, processed = [], [], []

    for i, at in enumerate(tqdm(frames, desc="Running xTB")):
        try:
            at.calc = xtb_calc
            e = at.get_potential_energy()
            f = at.get_forces()

            at.info["energy_xtb"] = e
            at.arrays["forces_xtb"] = f

            energies.append(e)
            forces.append(f)
            processed.append(at)
            logging.info(f"✅ Frame {i} | Energy = {e:.6f} eV")
        except Exception as e:
            logging.error(f"❌ Frame {i} FAILED | {e}")
            continue

    # 결과 저장
    if args.out_xyz:
        write(args.out_xyz, processed)
        logging.info(f"✅ Saved result XYZ → {args.out_xyz}")

    if args.out_npz:
        np.savez(args.out_npz, energy=np.array(energies), forces=np.array(forces))
        logging.info(f"✅ Saved NPZ dataset → {args.out_npz}")

    print("\n🎉 xTB Evaluation Complete")
    print(f"   ✅ Frames processed: {len(processed)}")
    print(f"   📄 Output XYZ: {args.out_xyz}")
    print(f"   💾 Output NPZ: {args.out_npz}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="xTB energy/force evaluator using ASE interface.")

    parser.add_argument("--input", type=str, default="xtbopt/all_ionpairs.xyz",
                        help="Input multi-frame XYZ file")
    parser.add_argument("--out_xyz", type=str, default="xtbopt/xtb_results.xyz",
                        help="Output XYZ file with xtb energy/forces")
    parser.add_argument("--out_npz", type=str, default="xtbopt/xtb_results.npz",
                        help="Output NPZ dataset file")
    parser.add_argument("--method", type=str, default="GFN2-xTB",
                        help="xTB method (GFN2-xTB / GFN1-xTB / GFN-FF)")
    parser.add_argument("--log", type=str, default="xtb_evaluate.log",
                        help="Log file name")

    args = parser.parse_args()
    setup_logging(args.log)
    run_xtb(args)

