import os
import argparse
import logging
from tqdm import tqdm
from ase.io import read, write

def collect_xyz_files(packmol_dir):
    xyz_files = []
    for case in sorted(os.listdir(packmol_dir)):
        case_dir = os.path.join(packmol_dir, case)
        xyz_path = os.path.join(case_dir, "ionpairs.xyz")
        if os.path.isfile(xyz_path):
            xyz_files.append(xyz_path)
        else:
            logging.warning(f"⚠️ ionpairs.xyz 없음: {case_dir}")
    return xyz_files

def main(args):
    packmol_dir = args.packmol_dir
    output_xyz = args.output

    logging.info("✅ XYZ 병합 시작")
    xyz_files = collect_xyz_files(packmol_dir)

    if not xyz_files:
        logging.error("❌ XYZ 파일이 없습니다. 스크립트 종료.")
        exit()

    all_frames = []
    for xyz_file in tqdm(xyz_files, desc="Merging xyz files"):
        try:
            frames = read(xyz_file, ":")  # 모든 frame 읽기
            all_frames.extend(frames)
            logging.info(f"✅ 추가됨: {xyz_file} ({len(frames)} 프레임)")
        except Exception as e:
            logging.error(f"❌ 읽기 실패: {xyz_file} | {e}")

    write(output_xyz, all_frames)
    logging.info(f"✅ 병합 완료 → {output_xyz}")
    print(f"✅ 병합된 XYZ 저장 완료: {output_xyz} (총 {len(all_frames)} 프레임)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge ionpairs.xyz files into one multi-frame XYZ file.")
    parser.add_argument(
        "--packmol_dir",
        type=str,
        default="packmol/packmol_cases",
        help="Packmol case 디렉토리 경로 (기본값: packmol/packmol_cases)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="all_ionpairs.xyz",
        help="출력 멀티프레임 XYZ 파일 이름 (기본값: all_ionpairs.xyz)"
    )
    parser.add_argument(
        "--log",
        type=str,
        default="merge_xyz.log",
        help="로그 파일 이름 (기본값: merge_xyz.log)"
    )
    args = parser.parse_args()

    # 로그 설정
    logging.basicConfig(
        filename=args.log,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    main(args)

