python 00_md.py \
  --data_path /data/dydtkddhkdwk/20251022_MACE_MD_LAMMPS/00_bulkmake/packmol_cases/case_001/ionpairs.xyz \
  --model_path /data/dydtkddhkdwk/20251021_MACE_FIT/20251021_try03_cuequivariance_test_with_Foundation_MACEOFF/MACE_OFF_models/mace_off_finetune_stagetwo_compiled.model \
  --output_md /data/dydtkddhkdwk/20251022_MACE_MD_LAMMPS/01_LAMMPS_MD/00_md_gpu/mace01_md.xyz \
  --temp 1200 \
  --steps 2000 \
  --interval 10 \
  --device cuda \
  --log_file /data/dydtkddhkdwk/20251022_MACE_MD_LAMMPS/01_LAMMPS_MD/00_md_gpu/md_run_gpu.log
