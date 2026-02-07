import warnings
warnings.filterwarnings("ignore")
from mace.cli.run_train import main as mace_train_main
import sys
import logging
import argparse
import os

def expand_config(base_path, out_path):
    """
    base.yml을 읽어서 ${WORKDIR}를 절대경로로 치환 후 out_path로 저장
    """
    workdir = os.getcwd()
    with open(base_path, "r") as f:
        text = f.read()

    # ${WORKDIR} 또는 $WORKDIR 치환
    text = text.replace("${WORKDIR}", workdir).replace("$WORKDIR", workdir)

    with open(out_path, "w") as f:
        f.write(text)

    print(f"✅ Config expanded and saved to {out_path}")
    return out_path

def train_mace(config_file):
    logging.getLogger().handlers.clear()
    sys.argv = ["mace-train", "--config", config_file]
    mace_train_main()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=str, required=True, help="Path to base YAML config file (with WORKDIR placeholders)")
    parser.add_argument("--config", type=str, default="config.yml", help="Output full config file path")
    args = parser.parse_args()

    print(f"🔥 Loading base config: {args.base}")
    config_path = expand_config(args.base, args.config)

    print(f"🚀 Starting MACE training with {config_path}")
    train_mace(config_path)
    print(f"✅ Training completed successfully.")

