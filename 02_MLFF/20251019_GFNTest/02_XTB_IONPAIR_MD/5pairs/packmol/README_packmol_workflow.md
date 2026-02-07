
# VBBI+/TFSI- Ion-Pair Packmol Workflow

This workflow builds many Packmol cases that place **5 TFSI** and **5 VBBI** conformers
(per case) randomly inside a box. It **samples conformers randomly from multi-frame XYZ files**.

## Files
- `build_packmol_ionpairs.py` — main script (logging + tqdm).
- `template_base_single.inp` — example Packmol template (manual use).
- This README.

## Requirements
- Python 3.8+
- `tqdm`
- Packmol executable in PATH (or pass `--packmol /path/to/packmol`).

## Usage (example)

```bash
python build_packmol_ionpairs.py   --tfsifiles ../../TFSI_confcross.xyz ../../TFSI_crest_conformers.xyz ../../TFSI_crest_rotamers.xyz   --vbbifiles ../../VBBI_confcross.xyz ../../VBBI_crest_conformers.xyz ../../VBBI_crest_rotamers.xyz   --ncases 50   --ntfsi 5 --nvbbi 5   --box 0 0 0 80 80 80   --tolerance 2.0   --outroot ./packmol_cases   --packmol packmol   --seed 1000   --timeout 120   --verbose
```

- Outputs will be inside `./packmol_cases/case_XXX/ionpairs.xyz`.
- Each case uses a distinct RNG seed (`seed + i`).

## Tips
- If Packmol often fails, enlarge the box (e.g., `--box 0 0 0 100 100 100`) or increase `--tolerance`.
- Keep `--keep_temps` to retain per-case `structures/` for debugging.
- Use `--dry_run` to only create inputs without running Packmol.
