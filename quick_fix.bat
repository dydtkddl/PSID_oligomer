@echo off
echo ========================================
echo Git Large Files - Quick Clean Solution
echo ========================================

echo [1/5] Creating backup branch...
git branch backup_before_cleanup

echo [2/5] Removing .git folder...
rmdir /s /q .git

echo [3/5] Reinitializing Git...
git init

echo [4/5] Setting up .gitignore...
(
echo # === Large simulation files ===
echo 20251208_MS/Oligomer_COMPASS_Files/Documents/
echo *.trj
echo *.mol2
echo *.xtd
echo *.arc
echo *.xtc
echo *.dcd
echo *.lammpstrj
echo.
echo # === Backup folders ===
echo *_BACKUP*/
) > .gitignore

echo [5/5] Adding files and pushing...
git add .
git commit -m "Clean repository: Code and analysis scripts only

- Excluded all large simulation files (*.trj, *.mol2, etc.)
- Focused on code, scripts, and documentation
- Large files preserved locally in Documents folder"

git branch -M main
git remote add origin https://github.com/dydtkddl/PSID_oligomer.git
git push -u origin main --force

echo.
echo ========================================
echo SUCCESS! Repository cleaned and pushed
echo ========================================
echo Large files are still available locally in:
echo %CD%\20251208_MS\Oligomer_COMPASS_Files\Documents\
echo.
pause