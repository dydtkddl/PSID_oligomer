import warnings
warnings.filterwarnings("ignore")
from mace.cli.run_train import main as mace_train_main
import sys
import logging
import argparse

def train_mace(config_file):
    logging.getLogger().handlers.clear()
    sys.argv = ["mace-train", "--config", config_file]
    mace_train_main()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, required=True, help="Path to YAML config file")
    args = parser.parse_args()

    print(f"🔥 Starting MACE training with config: {args.config}")
    train_mace(args.config)

