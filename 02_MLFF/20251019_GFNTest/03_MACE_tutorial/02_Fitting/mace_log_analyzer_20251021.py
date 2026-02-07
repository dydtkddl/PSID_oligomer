#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mace_log_analyzer.py

- MACE 학습 로그에서 'TRAINING' 섹션만 정확히 추출
- 에폭별 loss, RMSE_E_per_atom, RMSE_F, 타임스탬프 파싱
- 초기(Initial) 지표도 별도 컬럼로 저장
- 75 epoch 기준으로 스테이지 전/후 영역 시각적 구분 (음영/세로선)
- 에폭별 소요시간, 누적시간 계산
- CSV/PNG 결과 저장
- argparse, tqdm, logging 사용

사용 예)
python mace_log_analyzer.py --log logs/train_99443.log --outdir results --stage2_epoch 75
"""

import re
import os
import sys
import csv
import math
import argparse
import logging
from datetime import datetime
from typing import List, Dict, Any, Tuple, Optional

from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt


# --------- 정규식 패턴 ---------
TS_PATTERN = r"(?P<ts>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3})"
INITIAL_PATTERN = re.compile(
    rf"{TS_PATTERN} .*?Initial:.*?loss=(?P<loss>[-+eE0-9\.]+),\s*RMSE_E_per_atom=\s*(?P<rmse_e>[-+eE0-9\.]+)\s*meV,\s*RMSE_F=\s*(?P<rmse_f>[-+eE0-9\.]+)\s*meV\s*/\s*A"
)
EPOCH_PATTERN = re.compile(
    rf"{TS_PATTERN} .*?Epoch\s+(?P<epoch>\d+):.*?loss=(?P<loss>[-+eE0-9\.]+),\s*RMSE_E_per_atom=\s*(?P<rmse_e>[-+eE0-9\.]+)\s*meV,\s*RMSE_F=\s*(?P<rmse_f>[-+eE0-9\.]+)\s*meV\s*/\s*A"
)
TRAIN_BEGIN_PATTERN = re.compile(r"=+TRAINING=+")
TRAIN_STARTED_PATTERN = re.compile(r"Started training")
TRAIN_END_PATTERN = re.compile(r"Training complete")

# --------- 유틸 ---------
def parse_ts(s: str) -> datetime:
    # 예: "2025-10-21 15:50:46.439"
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S.%f")

def find_training_sections(lines: List[str]) -> List[Tuple[int, int]]:
    """
    '===========TRAINING===========' ~ 'Training complete' 구간 인덱스 목록 반환
    """
    sections = []
    i = 0
    n = len(lines)
    while i < n:
        # TRAINING 헤더 찾기
        while i < n and not TRAIN_BEGIN_PATTERN.search(lines[i]):
            i += 1
        if i >= n:
            break
        start_idx = i
        # 'Training complete'까지
        while i < n and not TRAIN_END_PATTERN.search(lines[i]):
            i += 1
        if i < n:
            end_idx = i  # 'Training complete' 라인 포함
            sections.append((start_idx, end_idx))
        i += 1
    return sections

def parse_training_block(lines: List[str]) -> Tuple[pd.DataFrame, Optional[Dict[str, Any]]]:
    """
    TRAINING 블록에서 Initial, Epoch 라인 파싱 → DataFrame 반환
    """
    rows = []
    initial_row = None

    for line in tqdm(lines, desc="Parsing training block", unit="line"):
        m0 = INITIAL_PATTERN.search(line)
        if m0:
            ts = parse_ts(m0.group("ts"))
            initial_row = {
                "is_initial": True,
                "epoch": -1,
                "timestamp": ts,
                "loss": float(m0.group("loss")),
                "rmse_e_mev_per_atom": float(m0.group("rmse_e")),
                "rmse_f_mev_per_A": float(m0.group("rmse_f")),
            }
            # Initial도 rows에 넣어두면 누적 시간 계산이 편함(에폭 -1로 표시)
            rows.append(initial_row)
            continue

        m = EPOCH_PATTERN.search(line)
        if m:
            ts = parse_ts(m.group("ts"))
            rows.append({
                "is_initial": False,
                "epoch": int(m.group("epoch")),
                "timestamp": ts,
                "loss": float(m.group("loss")),
                "rmse_e_mev_per_atom": float(m.group("rmse_e")),
                "rmse_f_mev_per_A": float(m.group("rmse_f")),
            })

    if not rows:
        raise ValueError("훈련 블록에서 유효한 Initial/Epoch 라인을 찾지 못했습니다.")

    df = pd.DataFrame(rows).sort_values(["epoch", "timestamp"]).reset_index(drop=True)

    # 에폭별 소요시간 및 누적시간 계산
    # Initial(-1) → Epoch 0 → 1 → ...
    df["epoch_duration_sec"] = None
    df["cumulative_sec"] = None

    prev_ts = None
    cumulative = 0.0
    for idx, r in df.iterrows():
        ts = r["timestamp"]
        if prev_ts is None:
            dur = 0.0
        else:
            dur = (ts - prev_ts).total_seconds()
        cumulative += dur
        df.at[idx, "epoch_duration_sec"] = dur
        df.at[idx, "cumulative_sec"] = cumulative
        prev_ts = ts

    return df, initial_row

def plot_metrics(df: pd.DataFrame, out_png: str, stage2_epoch: int):
    """
    Loss, RMSE_E_per_atom, RMSE_F vs Epoch (Initial 제외)
    Stage2 경계(세로선) 및 영역 음영 표시
    """
    data = df[df["epoch"] >= 0].copy()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(data["epoch"], data["loss"], label="Loss")
    ax.plot(data["epoch"], data["rmse_e_mev_per_atom"], label="RMSE_E_per_atom (meV/atom)")
    ax.plot(data["epoch"], data["rmse_f_mev_per_A"], label="RMSE_F (meV/Å)")

    # Stage shading
    min_ep = int(data["epoch"].min()) if len(data) else 0
    max_ep = int(data["epoch"].max()) if len(data) else stage2_epoch
    ax.axvspan(min_ep - 0.5, stage2_epoch - 0.5, alpha=0.12)
    ax.axvspan(stage2_epoch - 0.5, max_ep + 0.5, alpha=0.08)
    ax.axvline(stage2_epoch - 0.5, linestyle="--")

    ax.set_xlabel("Epoch")
    ax.set_ylabel("Metric value")
    ax.set_title("Validation metrics over epochs")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

def plot_time(df: pd.DataFrame, out_png: str, stage2_epoch: int):
    """
    Cumulative train time vs Epoch (Initial 제외)
    """
    data = df[df["epoch"] >= 0].copy()
    # 동일 epoch가 하나만 있다고 가정(표준 로그). 혹시 중복 있으면 마지막 ts 사용.
    data = data.drop_duplicates(subset=["epoch"], keep="last").sort_values("epoch")

    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.plot(data["epoch"], data["cumulative_sec"], marker="o", linewidth=1.5)
    # Stage shading
    min_ep = int(data["epoch"].min()) if len(data) else 0
    max_ep = int(data["epoch"].max()) if len(data) else stage2_epoch
    ax.axvspan(min_ep - 0.5, stage2_epoch - 0.5, alpha=0.12)
    ax.axvspan(stage2_epoch - 0.5, max_ep + 0.5, alpha=0.08)
    ax.axvline(stage2_epoch - 0.5, linestyle="--")

    ax.set_xlabel("Epoch")
    ax.set_ylabel("Cumulative time (sec)")
    ax.set_title("Cumulative training time")
    ax.grid(True, linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

def save_csv(df: pd.DataFrame, out_csv: str):
    # timestamp를 문자열로
    df2 = df.copy()
    df2["timestamp"] = df2["timestamp"].dt.strftime("%Y-%m-%d %H:%M:%S.%f").str[:-3]
    df2.to_csv(out_csv, index=False)

def main():
    ap = argparse.ArgumentParser(description="MACE training log analyzer (TRAINING block only)")
    ap.add_argument("--log", required=True, help="MACE training log file path")
    ap.add_argument("--outdir", default="mace_log_report", help="Output directory")
    ap.add_argument("--stage2_epoch", type=int, default=75, help="Stage 2 starts at this epoch")
    ap.add_argument("--section_index", type=int, default=-1,
                    help="Which TRAINING section to parse (default: last, -1)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ---- logging ----
    log_path = os.path.join(args.outdir, "analyzer.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode="w", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logging.info("Start MACE log analysis")
    logging.info(f"log file: {args.log}")
    logging.info(f"outdir:   {args.outdir}")

    # ---- read file ----
    with open(args.log, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    sections = find_training_sections(lines)
    if not sections:
        logging.error("TRAINING 섹션을 찾지 못했습니다.")
        sys.exit(1)

    idx = args.section_index if args.section_index >= 0 else (len(sections) - 1)
    if idx < 0 or idx >= len(sections):
        logging.error(f"section_index 범위 오류: {idx}, 섹션 수: {len(sections)}")
        sys.exit(1)

    s, e = sections[idx]
    sel_lines = lines[s:e+1]
    logging.info(f"선택된 TRAINING 섹션 인덱스: {idx}, 라인 범위: [{s}, {e}] (총 {len(sel_lines)}줄)")

    # ---- parse ----
    df, initial_row = parse_training_block(sel_lines)
    logging.info(f"파싱 완료: 총 {len(df)} 레코드 (Initial 포함). 에폭 범위: {df['epoch'].min()} ~ {df['epoch'].max()}")

    # ---- save CSV ----
    out_csv = os.path.join(args.outdir, "metrics_epoch_times.csv")
    save_csv(df, out_csv)
    logging.info(f"CSV 저장: {out_csv}")

    # ---- plots ----
    out_png_metrics = os.path.join(args.outdir, "plot_metrics.png")
    out_png_time = os.path.join(args.outdir, "plot_time.png")

    plot_metrics(df, out_png_metrics, args.stage2_epoch)
    logging.info(f"메트릭 플롯 저장: {out_png_metrics}")

    plot_time(df, out_png_time, args.stage2_epoch)
    logging.info(f"시간 플롯 저장: {out_png_time}")

    # ---- 요약 ----
    # 최종(가장 큰 epoch) 지표 로그
    last_ep_rows = df[df["epoch"] == df["epoch"].max()]
    if not last_ep_rows.empty:
        r = last_ep_rows.iloc[-1]
        logging.info(
            "최종 에폭 요약: epoch=%s, loss=%.6f, RMSE_E=%.3f meV/atom, RMSE_F=%.3f meV/Å, 누적시간=%.1f sec",
            int(r["epoch"]), r["loss"], r["rmse_e_mev_per_atom"], r["rmse_f_mev_per_A"], r["cumulative_sec"]
        )

    logging.info("완료")


if __name__ == "__main__":
    main()

