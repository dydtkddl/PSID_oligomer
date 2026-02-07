
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build Packmol inputs for VBBI+/TFSI- ion-pair boxes by sampling random conformers.

Features
- Reads multi-frame XYZ files for TFSI and VBBI (you can pass multiple files per species)
- Randomly samples ntfsi and nvbbi conformers per case (default 5 and 5)
- Generates per-case Packmol input with one 'structure' block per selected frame (number 1 each)
- Runs Packmol with a timeout and logs everything
- Saves outputs into structured folders (one directory per case)
- Provides tqdm progress bar and rich logging (console + file)

Example
-------
python build_packmol_ionpairs.py \  --tfsifiles ../../TFSI_confcross.xyz ../../TFSI_crest_conformers.xyz ../../TFSI_crest_rotamers.xyz \  --vbbifiles ../../VBBI_confcross.xyz ../../VBBI_crest_conformers.xyz ../../VBBI_crest_rotamers.xyz \  --ncases 20 \  --ntfsi 5 --nvbbi 5 \  --box 0 0 0 80 80 80 \  --tolerance 2.0 \  --outroot ./packmol_cases \  --packmol packmol \  --seed 1234 \  --timeout 120 \  --verbose

Notes
-----
- Packmol places structures with the given tolerance; if it fails often, try enlarging the box or tolerance.
- For different conformers per molecule instance, we generate one structure block per selected frame.
- If your multi-frame files are huge, consider pre-filtering (e.g., top-N by energy) to speed up sampling.
"""

import os
import sys
import random
import shutil
import logging
import argparse
import subprocess
from pathlib import Path
from typing import List, Tuple
from tqdm import tqdm

# ---------------
# Utilities
# ---------------

def setup_logger(log_path: Path, verbose: bool) -> logging.Logger:
    logger = logging.getLogger("packmol_builder")
    logger.setLevel(logging.DEBUG)
    # Clear handlers if rerun in same interpreter
    if logger.handlers:
        for h in list(logger.handlers):
            logger.removeHandler(h)

    fmt = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    # File handler
    fh = logging.FileHandler(log_path, mode='w', encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    return logger
def read_multiframe_xyz(xyz_path: Path) -> List[str]:
    """
    Robust XYZ reader: safely parses multi-frame XYZ files.
    Handles missing newline, extra spaces, and avoids trimming last atom.
    """
    frames = []
    with xyz_path.open('r', encoding='utf-8', errors='ignore') as f:
        lines = [l.rstrip("\n") for l in f.readlines()]

    i = 0
    nlines = len(lines)

    while i < nlines:
        # Skip empty lines
        while i < nlines and not lines[i].strip():
            i += 1
        if i >= nlines:
            break

        # Try reading natoms
        try:
            natoms = int(lines[i].strip())
        except ValueError:
            i += 1
            continue

        # Check frame completeness
        if i + 1 + natoms > nlines:
            # Incomplete block, abort
            break

        # Extract block safely
        start = i
        comment_line = lines[i + 1]
        atom_lines = lines[i + 2:i + 2 + natoms]

        # Reconstruct XYZ block
        block = "\n".join(
            [str(natoms), comment_line] + atom_lines
        ) + "\n"
        frames.append(block)

        # Move index to next possible frame
        i = i + 2 + natoms

    return frames


def write_single_xyz(block: str, out_path: Path) -> None:
    out_path.write_text(block, encoding='utf-8')


def build_packmol_input(case_dir: Path,
                        tfsis: List[Path],
                        vbbis: List[Path],
                        box: Tuple[float, float, float, float, float, float],
                        tolerance: float,
                        seed: int,
                        output_name: str = "ionpairs.xyz") -> Path:
    """
    Create a Packmol input file with 10 structure blocks (5 TFSI + 5 VBBI by default).
    Each 'structure' references a single-frame XYZ file and 'number 1'.
    """
    x1,y1,z1,x2,y2,z2 = box
    lines = []
    lines.append(f"tolerance {tolerance}")
    lines.append("filetype xyz")
    lines.append(f"output {output_name}")
    lines.append("")

    # TFSI blocks
    for p in tfsis:
        lines.append(f"structure {p.as_posix()}")
        lines.append("  number 1")
        lines.append(f"  inside box {x1} {y1} {z1} {x2} {y2} {z2}")
        lines.append("end structure")
        lines.append("")

    # VBBI blocks
    for p in vbbis:
        lines.append(f"structure {p.as_posix()}")
        lines.append("  number 1")
        lines.append(f"  inside box {x1} {y1} {z1} {x2} {y2} {z2}")
        lines.append("end structure")
        lines.append("")

    # seed & random placement
    lines.append(f"seed {seed}")
    lines.append("randominitialpoint")
    content = "\n".join(lines)

    inp_path = case_dir / "packmol.inp"
    inp_path.write_text(content, encoding='utf-8')
    return inp_path


def run_packmol(packmol_exe: str, inp_path: Path, work_dir: Path,
                timeout_s: int, logger: logging.Logger) -> int:
    cmd = f"{packmol_exe} < {inp_path.name}"
    logger.debug(f"Running in {work_dir}: {cmd}")
    try:
        result = subprocess.run(cmd, cwd=work_dir, shell=True, timeout=timeout_s)
        return result.returncode
    except subprocess.TimeoutExpired:
        logger.error(f"Packmol timeout after {timeout_s} s")
        return 124  # common timeout code


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


# ---------------
# Main
# ---------------

def main():
    ap = argparse.ArgumentParser(description="Build Packmol inputs for VBBI+/TFSI- ion pairs by sampling conformers.")
    ap.add_argument("--tfsifiles", nargs="+", required=True,
                    help="Paths to multi-frame XYZ files for TFSI (one or more).")
    ap.add_argument("--vbbifiles", nargs="+", required=True,
                    help="Paths to multi-frame XYZ files for VBBI (one or more).")
    ap.add_argument("--ncases", type=int, default=10, help="Number of Packmol cases to generate.")
    ap.add_argument("--ntfsi", type=int, default=5, help="Number of TFSI molecules per case.")
    ap.add_argument("--nvbbi", type=int, default=5, help="Number of VBBI molecules per case.")
    ap.add_argument("--box", nargs=6, type=float, default=[0,0,0,80,80,80],
                    help="Box as 6 floats: x1 y1 z1 x2 y2 z2 (default 0 0 0 80 80 80).")
    ap.add_argument("--tolerance", type=float, default=2.0, help="Packmol tolerance (Å).")
    ap.add_argument("--outroot", type=str, default="packmol_cases", help="Root directory for generated cases.")
    ap.add_argument("--packmol", type=str, default="packmol", help="Packmol executable path.")
    ap.add_argument("--timeout", type=int, default=120, help="Timeout for packmol execution (seconds).")
    ap.add_argument("--seed", type=int, default=1234, help="Base RNG seed. Each case uses seed+i.")
    ap.add_argument("--keep_temps", action="store_true", help="Keep per-case extracted frame files.")
    ap.add_argument("--dry_run", action="store_true", help="Only create inputs, do not run packmol.")
    ap.add_argument("--verbose", action="store_true", help="Verbose console logging.")
    args = ap.parse_args()

    outroot = Path(args.outroot).resolve()
    ensure_dir(outroot)
    log_path = outroot / "build.log"
    logger = setup_logger(log_path, args.verbose)

    # Read frames
    tfsifiles = [Path(p).resolve() for p in args.tfsifiles]
    vbbifiles = [Path(p).resolve() for p in args.vbbifiles]

    logger.info("Reading multi-frame XYZ files...")
    tfsiframes = []
    for p in tfsifiles:
        frames = read_multiframe_xyz(p)
        logger.info(f"TFSI file {p.name}: {len(frames)} frames")
        tfsiframes.extend(frames)
    vbbiframes = []
    for p in vbbifiles:
        frames = read_multiframe_xyz(p)
        logger.info(f"VBBI file {p.name}: {len(frames)} frames")
        vbbiframes.extend(frames)

    if len(tfsiframes) == 0 or len(vbbiframes) == 0:
        logger.error("No frames read. Check file paths and formats.")
        sys.exit(2)

    logger.info(f"Total TFSI frames: {len(tfsiframes)}; total VBBI frames: {len(vbbiframes)}")

    # RNG
    rng = random.Random(args.seed)

    # Process cases
    success, failed, timeouts = 0, 0, 0
    for i in tqdm(range(1, args.ncases + 1), desc="Cases"):
        case_seed = args.seed + i
        case_dir = outroot / f"case_{i:03d}"
        struct_dir = case_dir / "structures"
        ensure_dir(case_dir)
        ensure_dir(struct_dir)

        # Sample frames (without replacement if enough, else with replacement)
        def sample_blocks(pool, k):
            if len(pool) >= k:
                idxs = rng.sample(range(len(pool)), k)
            else:
                idxs = [rng.randrange(len(pool)) for _ in range(k)]
            return idxs

        t_idx = sample_blocks(tfsiframes, args.ntfsi)
        v_idx = sample_blocks(vbbiframes, args.nvbbi)
        logger.debug(f"Case {i:03d}: TFSI idx {t_idx}; VBBI idx {v_idx}")

        # Write per-frame XYZ files
        t_paths = []
        for j, idx in enumerate(t_idx, start=1):
            p = struct_dir / f"TFSI_{j:02d}.xyz"
            write_single_xyz(tfsiframes[idx], p)
            t_paths.append(p)
        v_paths = []
        for j, idx in enumerate(v_idx, start=1):
            p = struct_dir / f"VBBI_{j:02d}.xyz"
            write_single_xyz(vbbiframes[idx], p)
            v_paths.append(p)

        # Build packmol input
        inp_path = build_packmol_input(case_dir, t_paths, v_paths,
                                       tuple(args.box), args.tolerance,
                                       case_seed, output_name="ionpairs.xyz")

        # Run packmol
        if not args.dry_run:
            rc = run_packmol(args.packmol, inp_path, case_dir, args.timeout, logger)
            if rc == 0:
                success += 1
                logger.info(f"[OK] case_{i:03d} completed.")
            #    if not args.keep_temps:
                    # Clean structure files to save space
    #                try:
     #                   shutil.rmtree(struct_dir)
      #              except Exception as e:
       #                 logger.warning(f"Could not remove temp dir: {struct_dir} ({e})")
            elif rc == 124:
                timeouts += 1
                logger.error(f"[TIMEOUT] case_{i:03d}")
            else:
                failed += 1
                logger.error(f"[FAIL] case_{i:03d} (rc={rc})")
        else:
            logger.info(f"[DRY] case_{i:03d} input prepared (seed={case_seed}).")

    logger.info("=== Summary ===")
    logger.info(f"Success: {success} / {args.ncases}")
    logger.info(f"Timeouts: {timeouts}")
    logger.info(f"Failed:   {failed}")
    print(f"\nSummary -> Success: {success}, Timeouts: {timeouts}, Failed: {failed}")
    print(f"Logs: {log_path}")

if __name__ == "__main__":
    main()
