import os
import math
import argparse
import pandas as pd
from joblib import Parallel, delayed, cpu_count


# ✅ row를 GPT-friendly 텍스트로 변환
def format_row_for_prompt(row):
    return f"""
Atom Type Data:
- Atom index: {row['Atom_Index']}
- Force field raw info: {row['Info']}
- Partial charge: {row['Charge']}
- Atomic mass: {row['DataMass']}
- Element guess: {row['ElementGuess']}
- Atomic number: {row['AtomicNumber']}
- Source file: {row['Source']}
"""


# ✅ System Prompt
SYSTEM_PROMPT = """
You are a cheminformatics expert skilled in force field atom typing, SMILES/SMARTS pattern matching, and chemical environment classification.
You will convert textual OPLS-AA atom-type descriptions into structured chemistry metadata.

Your task:
- Interpret the chemical environment described.
- Identify the exact atom location described (e.g. “C2 in CH2CN” means the carbon alpha to the nitrile group).
- Distinguish between two SMARTS patterns:
  1. A general SMARTS describing the functional group environment.
  2. A specific SMARTS that matches only the exact atom described (atomic-focused).
- Infer hybridization and nearest neighbor atom types when possible.

Output Format (JSON):
Return a valid JSON object with these keys:
- atom_type (integer)
- element (string)
- chemical_environment (string)
- functional_group (string)
- hybridization (string)
- example_molecule_smiles (string)
- group_smarts (string)
- atom_specific_smarts (string)
- description (string)
- uncertainty (string, 0~100 scale)

Rules:
- Do NOT include any explanation or commentary outside the JSON.
- SMARTS must be chemically meaningful and precise.
- If the atom index is indicated (e.g. C2), reflect that in SMARTS.
- If unsure, include "uncertainty": "medium".
"""


# ✅ Batch 처리 함수 (병렬 실행)
def process_batch(batch_df, batch_id, save_dir):
    from openai import OpenAI  # ✅ 병렬 안전 위해 함수 내부에서 import
    client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])  # ✅ 환경변수 방식 사용 권장

    os.makedirs(save_dir, exist_ok=True)

    rows_text = "\n".join(format_row_for_prompt(row) for _, row in batch_df.iterrows())
    user_prompt = f"Convert the following OPLS-AA atom type rows:\n{rows_text}"

    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": user_prompt}
    ]

    try:
        response = client.chat.completions.create(
            model="gpt-5",
            messages=messages
        )
        output_text = response.choices[0].message.content
    except Exception as e:
        output_text = f"API_ERROR: {str(e)}"

    batch_df = batch_df.copy()
    batch_df["prompt"] = user_prompt
    batch_df["gpt_output"] = output_text

    out_path = os.path.join(save_dir, f"batch_{batch_id}.csv")
    batch_df.to_csv(out_path, index=False)
    return out_path


# ✅ 전체 배치 실행
def run_batches(df, batch_size, n_cpus, save_dir):
    total_batches = math.ceil(len(df) / batch_size)
    print(f"🚀 시작: 총 {len(df)}개 row → {total_batches}개 batch로 처리 (CPU {n_cpus}개 병렬 사용)")

    batches = [
        (df.iloc[i * batch_size:(i + 1) * batch_size], i)
        for i in range(total_batches)
    ]

    results = Parallel(n_jobs=n_cpus)(
        delayed(process_batch)(batch_df, batch_id, save_dir)
        for batch_df, batch_id in batches
    )
    print("✅ 전체 완료!")
    return results


# ✅ CLI 설정
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OPLS-AA → GPT 구조 변환 Batch Pipeline")

    parser.add_argument("--input", required=True, help="입력 CSV 파일 경로")
    parser.add_argument("--save_dir", default="batch_outputs", help="배치 저장 폴더")
    parser.add_argument("--batch_size", type=int, default=20, help="Batch당 Row 개수")
    parser.add_argument("--n_cpus", type=int, default=cpu_count(), help="사용할 CPU 개수")

    args = parser.parse_args()

    DF = pd.read_csv(args.input)
    run_batches(DF, batch_size=args.batch_size, n_cpus=args.n_cpus, save_dir=args.save_dir)
