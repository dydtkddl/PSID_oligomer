import os
import argparse
import logging
import subprocess
from ase.io import read, write
from tqdm import tqdm

def setup_logging(logfile):
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

def run_xtb_cli_opt(output_dir, threads):
    """
    CLI xTB 실행 (input.xyz 기준)
    """
    cmd = [
        "xtb", "input.xyz",
        "--opt", "loose",
        "--gfn", "2",
        "--parallel", str(threads)
    ]
    try:
        with open(os.path.join(output_dir, "xtb.out"), "w") as out, \
             open(os.path.join(output_dir, "xtb.err"), "w") as err:
            subprocess.run(cmd, cwd=output_dir, stdout=out, stderr=err, check=True)
        logging.info(f"✅ xTB opt success: {output_dir}")
    except subprocess.CalledProcessError:
        logging.error(f"❌ xTB opt FAILED: {output_dir}")

def main(args):
    # 입력 불러오기
    frames = read(args.input, ":")
    os.makedirs(args.output_dir, exist_ok=True)

    for i, frame in enumerate(tqdm(frames, desc=" Running xTB optimization ")):
        frame_dir = os.path.join(args.output_dir, f"frame_{i:04d}")
        os.makedirs(frame_dir, exist_ok=True)

        # 프레임 저장
        write(os.path.join(frame_dir, "input.xyz"), frame)

        # 실행
        run_xtb_cli_opt(frame_dir, args.threads)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run xTB CLI geometry optimization on frames")
    parser.add_argument("--input", default="xtbopt/all_ionpairs.xyz", help="Input multi-frame XYZ file")
    parser.add_argument("--output_dir", default="xtbopt/xtb_opt", help="Output directory")
    parser.add_argument("--threads", default=4, type=int, help="Threads for xTB parallelization")
    parser.add_argument("--log", default="xtb_cli_opt.log", help="Log file name")

    args = parser.parse_args()
    setup_logging(args.log)
    main(args)
