import os
import time
import warnings
warnings.filterwarnings("ignore")

import os


import time
import warnings
warnings.filterwarnings("ignore")

import torch
import argparse
from ase.io import read
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
import subprocess
import random

# ====================== LOGGER ======================
def log(msg, logfile):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    with open(logfile, "a", encoding="utf-8") as f:
        f.write(line + "\n")

# ====================== GPU MONITOR ======================
def log_gpu_usage(logfile):
    try:
        gpu = subprocess.check_output(
            "nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits",
            shell=True
        ).decode().strip()
        log(f"[GPU] Memory used: {gpu} MB", logfile)
    except:
        log("[GPU] Unable to read GPU usage", logfile)

# ====================== MD FUNCTION ======================
def simpleMD(init_conf, temp, calc, fname, interval, total_steps, logfile, device):
    # 1. Calculator 설정
    init_conf.set_calculator(calc)

    # 2. 초기 속도 분포 세팅
    random.seed(42)
    MaxwellBoltzmannDistribution(init_conf, temperature_K=temp)
    Stationary(init_conf)
    ZeroRotation(init_conf)

    # 3. Langevin Thermostat
    dyn = Langevin(init_conf, 1.0 * units.fs, temperature_K=temp, friction=0.1)
    log(f"[MD] Langevin thermostat initialized (Target Temp={temp}K, dt=1fs)", logfile)
import torch
import argparse
from ase.io import read
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
import subprocess
import random
import matplotlib.pyplot as plt

import seaborn as sns  # ✅ seaborn 스타일을 matplotlib에 등록
# plt.style.use("seaborn-darkgrid")
# plt.style.use("seaborn-v0_8-darkgrid")
# ====================== LOGGER ======================
def log(msg, logfile):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    with open(logfile, "a", encoding="utf-8") as f:
        f.write(line + "\n")

# ====================== GPU MONITOR ======================
def log_gpu_usage(logfile):
    try:
        gpu = subprocess.check_output(
            "nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits",
            shell=True
        ).decode().strip()
        log(f"[GPU] Memory used: {gpu} MB", logfile)
    except:
        log("[GPU] Unable to read GPU usage", logfile)

# ====================== MD FUNCTION ======================
def simpleMD(init_conf, temp, calc, fname, interval, total_steps, logfile, device):
    init_conf.set_calculator(calc)
    random.seed(42)
    MaxwellBoltzmannDistribution(init_conf, temperature_K=temp)
    Stationary(init_conf)
    ZeroRotation(init_conf)

    dyn = Langevin(init_conf, 1.0 * units.fs, temperature_K=temp, friction=0.1)
    log(f"[MD] Langevin thermostat initialized (T={temp}K, dt=1fs)", logfile)

    if os.path.exists(fname):
        os.remove(fname)

    time_fs, temperatures, energies = [], [], []

    def save_step():
        dyn.atoms.write(fname, append=True)
        time_fs.append(dyn.get_time() / units.fs)
        temperatures.append(dyn.atoms.get_temperature())
        energies.append(dyn.atoms.get_potential_energy() / len(dyn.atoms))
        if len(time_fs) % 10 == 0:
            log(f"[MD] Step {len(time_fs)} | E={energies[-1]:.4f} eV/atom | T={temperatures[-1]:.1f}K", logfile)
            if device == "cuda":
                log_gpu_usage(logfile)

    dyn.attach(save_step, interval=interval)

    log(f"[MD] Simulation started (total_steps={total_steps}, interval={interval})", logfile)
    start = time.time()
    dyn.run(total_steps)
    elapsed = (time.time() - start) / 60
    log(f"[MD] Simulation finished in {elapsed:.2f} min → saved: {fname}", logfile)

    return time_fs, energies, temperatures

# ====================== PLOTTING FUNCTION ======================
def plot_md(time_fs, energies, temperatures, output_path):
    plt.style.use("seaborn-darkgrid")
    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    
    # Energy Plot
    ax[0].plot(time_fs, energies, linewidth=1.8)
    ax[0].set_ylabel("Energy (eV/atom)")
    ax[0].set_title("MD Simulation Results (MACE)")

    # Temperature Plot
    ax[1].plot(time_fs, temperatures, color="red", linewidth=1.8)
    ax[1].set_ylabel("Temperature (K)")
    ax[1].set_xlabel("Time (fs)")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

# ====================== MAIN ======================
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--data_path", required=True)
    parser.add_argument("--model_path", required=True)
    parser.add_argument("--output_md", required=True)
    parser.add_argument("--plot_file", default="md_plot.png")
    parser.add_argument("--temp", type=float, default=1200)
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--interval", type=int, default=10)
    parser.add_argument("--device", choices=["cuda", "cpu"], default="cuda")
    parser.add_argument("--num_threads", type=int, default=None)
    parser.add_argument("--log_file", default="md_run.log")

    args = parser.parse_args()

    if os.path.exists(args.log_file):
        os.remove(args.log_file)

    log("=== INITIALIZING MD RUN ===", args.log_file)

    if args.device == "cuda" and torch.cuda.is_available():
        device = "cuda"
        log(f"[DEVICE] Using GPU: {torch.cuda.get_device_name(0)}", args.log_file)
    else:
        device = "cpu"
        log("[DEVICE] Using CPU mode", args.log_file)
        if args.num_threads:
            os.environ["OMP_NUM_THREADS"] = str(args.num_threads)
            os.environ["MKL_NUM_THREADS"] = str(args.num_threads)
            torch.set_num_threads(args.num_threads)
            log(f"[CPU] Threads = {args.num_threads}", args.log_file)

    log(f"[LOAD] Structure: {args.data_path}", args.log_file)
    init_conf = read(args.data_path)

    log(f"[MODEL] Loading MACE model: {args.model_path}", args.log_file)
    from mace.calculators import MACECalculator
    mace_calc = MACECalculator(model_paths=[args.model_path], device=device, default_dtype="float32")

    os.makedirs(os.path.dirname(args.output_md), exist_ok=True)
    time_fs, energies, temperatures = simpleMD(init_conf, args.temp, mace_calc,
                                               args.output_md, args.interval, args.steps,
                                               args.log_file, device)

    log(f"[PLOT] Saving MD plot to {args.plot_file}", args.log_file)
    plot_md(time_fs, energies, temperatures, args.plot_file)

    log("[DONE] All tasks complete ✅", args.log_file)

if __name__ == "__main__":
    main()


