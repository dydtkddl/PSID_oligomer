import argparse
from ase.io import read
from x3dase.visualize import view_x3d_n

def visualize_trajectory(traj_path, every, start, end):
    """
    X3D 기반 분자동역학 trajectory 시각화 함수
    """
    # ASE로 trajectory 읽기
    print(f"[INFO] Loading trajectory: {traj_path}")
    traj = read(traj_path, index=':')
    
    # 슬라이싱 옵션 적용
    traj = traj[start:end:every]
    print(f"[INFO] Loaded {len(traj)} frames for visualization")

    # X3D 시각화
    print("[INFO] Launching 3D viewer...")
    view_x3d_n(traj)  # <-- Jupyter 혹은 localhost에서 3D 상호작용 가능

def main():
    parser = argparse.ArgumentParser(description="3D Visualization for MACE MD Trajectory using x3dase")

    parser.add_argument(
        "--traj",
        required=True,
        type=str,
        help="Path to trajectory XYZ file (e.g. mace01_md.xyz)"
    )
    parser.add_argument(
        "--every",
        type=int,
        default=1,
        help="Frame stride for visualization (default: 1)"
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="Start frame index (default: 0)"
    )
    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help="End frame index (default: entire trajectory)"
    )

    args = parser.parse_args()
    visualize_trajectory(args.traj, args.every, args.start, args.end)

if __name__ == "__main__":
    main()

