import argparse
import re
from datetime import datetime
from tqdm import tqdm
import logging

# -----------------------
# 로그 설정
# -----------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

def parse_standard_time(line):
    """기본 로그 타임스탬프 예: 2025-10-23 23:56:18"""
    try:
        match = re.match(r"(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})", line)
        if match:
            return datetime.strptime(match.group(1), "%Y-%m-%d %H:%M:%S")
    except:
        return None
    return None

def parse_slurm_end_time(line):
    """슬럼 로그 종료 시간 포맷: Thu 23 Oct 2025 11:40:06 PM KST"""
    try:
        match = re.match(r"([A-Za-z]{3} \d{2} [A-Za-z]{3} \d{4} \d{2}:\d{2}:\d{2} (AM|PM) [A-Za-z]{3})", line)
        if match:
            return datetime.strptime(match.group(1), "%a %d %b %Y %I:%M:%S %p %Z")
    except:
        return None
    return None

def extract_times(log_file):
    start_time = None
    md_start_time = None
    md_end_time = None

    with open(log_file, "r") as f:
        lines = f.readlines()

    next_is_end_time = False
    for line in tqdm(lines, desc="Parsing log"):
        # 시작 시간 (맨 처음 timestamp)
        if not start_time:
            t = parse_standard_time(line)
            if t:
                start_time = t

        # MD 시작 시각
        if ("DEBUG: Running mixed MD" in line or "0.0%" in line) and not md_start_time:
            md_start_time = parse_standard_time(line)

        # Job 끝나는 라인 감지
        if "Job Coeted" in line:
            next_is_end_time = True
            continue

        # Job Coeted 다음 줄 → 종료 타임
        if next_is_end_time:
            slurm_t = parse_slurm_end_time(line.strip())
            if slurm_t:
                md_end_time = slurm_t
                next_is_end_time = False

    return start_time, md_start_time, md_end_time

def main():
    parser = argparse.ArgumentParser(description="MACE MD Runtime Analyzer")
    parser.add_argument("--log", required=True, help="Path to MACE MD log file")
    args = parser.parse_args()

    logging.info(f"Reading log file: {args.log}")
    start, md_start, md_end = extract_times(args.log)

    if not start:
        logging.error("시작 시간 감지 실패")
        return
    if not md_start:
        logging.error("MD 시작 시각 감지 실패 (0% 또는 DEBUG marker 없음)")
        return
    if not md_end:
        logging.error("슬럼 종료 시간 감지 실패")
        return

    warmup_time = md_start - start
    md_time = md_end - md_start

    print("\n========== MACE MD Runtime Summary ==========")
    print(f"[1] Warmup time (Init → MD start):     {warmup_time}")
    print(f"[2] Pure MD runtime (0% → 100%):       {md_time}")
    print(f"[3] Total runtime (Init → End):        {md_end - start}")
    print("=============================================\n")


if __name__ == "__main__":
    main()

