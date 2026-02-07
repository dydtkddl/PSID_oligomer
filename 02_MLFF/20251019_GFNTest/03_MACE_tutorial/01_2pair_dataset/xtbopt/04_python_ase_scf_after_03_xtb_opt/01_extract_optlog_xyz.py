#!/usr/bin/env python3
import os
import re
import argparse
import logging
from typing import List
import numpy as np
from tqdm import tqdm

from ase import Atoms
from ase.io import write

# ------------------------- Logging -------------------------
def setup_logging(logfile: str):
    logging.basicConfig(
        filename=logfile,
        level=logging.INFO,
        filemode="w",
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    console = logging.StreamHandler()
    console.setLevel(logging.WARNING)
    console.setFormatter(logging.Formatter("%(levelname)s - %(message)s"))
    logging.getLogger("").addHandler(console)

# ------------------------- Parser -------------------------
def parse_xtb_xyz_log(path: str) -> List[Atoms]:
    """
    Parse xTB optimization log (xtbopt.log) containing repeated XYZ blocks:
      <n_atoms>
      energy: <val> gnorm: ...
      <sym> <x> <y> <z>
      ...
    Returns a list of ASE Atoms frames.
    """
    frames: List[Atoms] = []
    if not os.path.isfile(path):
        logging.warning(f"Missing file: {path}")
        return frames

    try:
        with open(path, "r") as f:
            lines = f.readlines()
    except Exception as e:
        logging.error(f"Failed to read {path}: {e}")
        return frames

    i = 0
    n_lines = len(lines)
    energy_re = re.compile(r"energy:\s*([-+]?[\d\.Ee+-]+)")

    while i < n_lines:
        # 1) number-of-atoms line
        line = lines[i].strip()
        i += 1
        if not line:
            continue
        try:
            natoms = int(line)
        except ValueError:
            # Not a frame start, skip
            continue

        if i >= n_lines:
            break

        # 2) comment line
        comment = lines[i].rstrip("\n")
        i += 1

        # 3) natoms coordinate lines
        if i + natoms > n_lines:
            logging.warning(f"Truncated block in {path} around line {i}")
            break

        symbols = []
        positions = np.zeros((natoms, 3), dtype=float)
        ok = True
        for j in range(natoms):
            parts = lines[i + j].split()
            if len(parts) < 4:
                ok = False
                break
            symbols.append(parts[0])
            try:
                positions[j, 0] = float(parts[1])
                positions[j, 1] = float(parts[2])
                positions[j, 2] = float(parts[3])
            except ValueError:
                ok = False
                break
        i += natoms

        if not ok:
            logging.warning(f"Malformed atom block in {path}, skipping a frame")
            continue

        at = Atoms(symbols=symbols, positions=positions)
        m = energy_re.search(comment)
        if m:
            try:
                at.info["energy_xtb"] = float(m.group(1))
            except Exception:
                pass
        at.info["comment"] = comment
        frames.append(at)

    return frames

# ------------------------- Sampling -------------------------
def sample_indices(n: int, interval: int) -> List[int]:
    """
    Return indices [0, interval, 2*interval, ...] and ALWAYS include last (n-1).
    """
    if n <= 0:
        return []
    step = max(1, interval)
    idx = list(range(0, n, step))
    if idx[-1] != n - 1:
        idx.append(n - 1)
    return idx

# ------------------------- Split Helper -------------------------
def split_frames(frames: List[Atoms], n_parts: int) -> List[List[Atoms]]:
    """
    Split frames into n_parts chunks as evenly as possible.
    """
    if n_parts <= 1 or len(frames) == 0:
        return [frames]
    # numpy.array_split spreads remainder nicely
    indices = np.array_split(np.arange(len(frames)), n_parts)
    return [[frames[i] for i in idx] for idx in indices if len(idx) > 0]

def ensure_parent_dir(path: str):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)

# ------------------------- Main -------------------------
def main(args):
    setup_logging(args.log)
    logging.info("Starting XTB trajectory sampling & merge")

    # Collect frame_* directories
    frame_dirs = [
        d for d in os.listdir(args.input_dir)
        if d.startswith("frame_") and os.path.isdir(os.path.join(args.input_dir, d))
    ]
    frame_dirs.sort()

    all_sampled: List[Atoms] = []

    for d in tqdm(frame_dirs, desc="Reading & sampling logs"):
        log_path = os.path.join(args.input_dir, d, "xtbopt.log")
        frames = parse_xtb_xyz_log(log_path)

        if not frames:
            logging.error(f"No frames parsed from {log_path}")
            continue

        idx = sample_indices(len(frames), args.interval)
        sampled = [frames[k] for k in idx]
        all_sampled.extend(sampled)

        logging.info(f"{d}: total={len(frames)}, sampled={len(sampled)} (interval={args.interval})")

        if args.save_each:
            per_out = os.path.join(args.input_dir, d, "xtb_sampled.xyz")
            try:
                write(per_out, sampled)
                logging.info(f"Saved per-folder sampled xyz → {per_out}")
            except Exception as e:
                logging.error(f"Failed to write {per_out}: {e}")

    # Nothing collected
    if not all_sampled:
        print("❌ No frames collected. Check logs and input paths.")
        logging.error("No frames collected. Nothing written.")
        return

    # Derive base for outputs
    base_out = args.output  # e.g., xtbopt/xtb_sampled.xyz or xtb_sampled.xyz
    root, ext = os.path.splitext(base_out)
    if ext.lower() not in (".xyz", ".extxyz"):
        # default to xyz if no/unknown extension
        root = base_out
        ext = ".xyz"

    # 1) Save ALL merged
    all_out = f"{root}_all{ext}"
    ensure_parent_dir(all_out)
    write(all_out, all_sampled)
    print(f"✅ Sampling complete. Total merged frames: {len(all_sampled)}")
    print(f"📄 ALL → {all_out}")
    logging.info(f"Merged sampled frames written to {all_out}")

    # 2) Optionally split into N parts
    if args.split and args.split > 0:
        chunks = split_frames(all_sampled, args.split)
        for i, chunk in enumerate(chunks, start=1):
            part_out = f"{root}_part_{i}{ext}"
            ensure_parent_dir(part_out)
            write(part_out, chunk)
            logging.info(f"Saved split part {i}/{len(chunks)} → {part_out}")
            print(f"   ▸ PART {i:02d} → {part_out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract XYZ frames from xTB xtbopt.log files, sample by interval, merge into one, and optionally split into N parts."
    )
    parser.add_argument(
        "--input_dir", type=str, default="xtbopt/03_xtb_opt",
        help="Directory containing frame_XXXX folders"
    )
    parser.add_argument(
        "--output", type=str, default="xtbopt/xtb_sampled.xyz",
        help="Base output path (will produce *_all.xyz and optionally *_part_k.xyz)"
    )
    parser.add_argument(
        "--interval", type=int, default=30,
        help="Sampling interval (last frame is always included)"
    )
    parser.add_argument(
        "--split", type=int, default=0,
        help="Number of parts to split the merged XYZ into (0 = no split)"
    )
    parser.add_argument(
        "--save_each", action="store_true",
        help="Also save per-folder sampled xyz as xtb_sampled.xyz"
    )
    parser.add_argument(
        "--log", type=str, default="xtb_sample.log",
        help="Log filename"
    )
    args = parser.parse_args()
    main(args)
