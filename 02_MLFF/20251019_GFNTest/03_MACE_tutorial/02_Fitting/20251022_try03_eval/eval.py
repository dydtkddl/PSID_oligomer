import os
import sys
import time
import warnings
warnings.filterwarnings("ignore")

import torch
from tqdm import tqdm
from mace.cli.eval_configs import main as mace_eval_configs_main

# ====================== Visualization 패키지 ======================
from ase.io import read
import matplotlib.pyplot as plt
from aseMolec import pltProps as pp
from aseMolec import extAtoms as ea
from aseMolec import anaAtoms as aa
import numpy as np

# ====================== GPU CHECK ======================
print("=====================================================")
if torch.cuda.is_available():
    print(f"[INFO] CUDA is available ✅ Using GPU: {torch.cuda.get_device_name(0)}")
    print(f"[INFO] CUDA Device Count: {torch.cuda.device_count()}")
else:
    print("[WARN] CUDA NOT available ⚠ Falling back to CPU")
print("=====================================================\n")

# ====================== SIMPLE LOGGER ======================
def log(msg):
    print(f"[INFO] {time.strftime('%Y-%m-%d %H:%M:%S')} | {msg}")

# ====================== OUTPUT DIR ======================
os.makedirs("tests/mace_eval", exist_ok=True)

# ====================== EVAL FUNCTION ======================
def eval_mace(configs, model, output, batch_size=10, device="cuda"):
    log(f"🚀 Evaluating dataset: {configs}")
    os.environ["CUDA_VISIBLE_DEVICES"] = "0"
    os.environ["MACE_DEFAULT_DEVICE"] = device

    for _ in tqdm(range(1), desc=f"Running eval on {os.path.basename(configs)}"):
        sys.argv = [
            "program",
            "--configs", configs,
            "--model", model,
            "--output", output,
            "--batch_size", str(batch_size),
            "--device", device
        ]
        try:
            mace_eval_configs_main()
        except Exception as e:
            log(f"❌ Evaluation error: {e}")
            raise e

    log(f"✅ Done: {configs}\n")

# ====================== RMSE Plot 함수 ======================
def plot_RMSEs(db, labs, title_tag):
    ea.rename_prop_tag(db, 'MACE_energy', 'energy_mace')
    ea.rename_prop_tag(db, 'MACE_forces', 'forces_mace')

    plt.figure(figsize=(9, 6), dpi=100)

    plt.subplot(1, 3, 1)
    pp.plot_prop(ea.get_prop(db, 'bind', '_xtb', True).flatten(),
                 ea.get_prop(db, 'bind', '_mace', True).flatten(),
                 title='Binding Energy (eV/atom)', labs=labs, rel=False)

    plt.subplot(1, 3, 2)
    pp.plot_prop(ea.get_prop(db, 'info', 'energy_xtb', True).flatten(),
                 ea.get_prop(db, 'info', 'energy_mace', True).flatten(),
                 title='Total Energy (eV/atom)', labs=labs, rel=False)

    plt.subplot(1, 3, 3)
    pp.plot_prop(np.concatenate(ea.get_prop(db, 'arrays', 'forces_xtb')).flatten(),
                 np.concatenate(ea.get_prop(db, 'arrays', 'forces_mace')).flatten(),
                 title='Forces (eV/Å)', labs=labs, rel=False)

    plt.suptitle(f"RMSE Comparison - {title_tag}")
    plt.tight_layout()
    save_path = f"tests/mace_eval/RMSE_{title_tag}.png"
    plt.savefig(save_path)
    log(f"📊 Saved RMSE plot: {save_path}")
    plt.close()

# ====================== 절대 경로 유지 ======================
BASE_DIR = "/data/dydtkddhkdwk/20251021_MACE_FIT"
TRAIN_CONFIG = f"{BASE_DIR}/20251021_try03_cuequivariance_test_with_Foundation_MACEOFF/dataset_train_0.5_.extxyz"
TEST_CONFIG  = f"{BASE_DIR}/20251021_try03_cuequivariance_test_with_Foundation_MACEOFF/dataset_test_0.5_.extxyz"
MODEL_PATH   = f"{BASE_DIR}/20251021_try03_cuequivariance_test_with_Foundation_MACEOFF/MACE_OFF_models/mace_off_finetune_stagetwo_compiled.model"

# ====================== RUN ======================
log("✅ START EVALUATION")
eval_mace(TRAIN_CONFIG, MODEL_PATH, "tests/mace_eval/train_eval.xyz")
eval_mace(TEST_CONFIG, MODEL_PATH, "tests/mace_eval/test_eval.xyz")

log("📥 Loading evaluation results...")
train_data = read("tests/mace_eval/train_eval.xyz", ":")
test_data  = train_data[:6] + read("tests/mace_eval/test_eval.xyz", ":")

# ====================== RMSE Visualization ======================
log("📊 Generating RMSE plots...")
plot_RMSEs(train_data, labs=["XTB", "MACE"], title_tag="TRAIN")
plot_RMSEs(test_data, labs=["XTB", "MACE"], title_tag="TEST")

# ====================== 추가: Trans/Rot/Vib 비교 ======================
log("⚙️ Extracting molecular properties for vibrational analysis...")
ea.rename_prop_tag(train_data, 'energy_xtb', 'energy')
ea.rename_prop_tag(train_data, 'forces_xtb', 'forces')
ea.rename_prop_tag(test_data, 'energy_xtb', 'energy')
ea.rename_prop_tag(test_data, 'forces_xtb', 'forces')
# db_xtb = train_data

# db_mace = read("tests/mace_eval/train_eval.xyz", ":")
db_xtb = test_data

db_mace = train_data[:6] + read("tests/mace_eval/test_eval.xyz", ":")
ea.rename_prop_tag(db_mace, 'MACE_energy', 'energy')
ea.rename_prop_tag(db_mace, 'MACE_forces', 'forces')

aa.extract_molecs(db_xtb, intra_inter=True)
aa.extract_molecs(db_mace, intra_inter=True)

plt.figure()
pp.plot_trans_rot_vib(db_xtb, db_mace, labs=['XTB', 'MACE'])
vib_save_path = "tests/mace_eval/VIB_Comparison.png"
plt.savefig(vib_save_path)
log(f"🎯 Saved Trans/Rot/Vib plot: {vib_save_path}")

log("🎉 ALL EVALUATION + VISUALIZATION COMPLETE ✅")
