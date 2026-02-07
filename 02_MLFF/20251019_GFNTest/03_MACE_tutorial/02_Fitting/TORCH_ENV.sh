# ✅ CUDA 12.1 대응 PyTorch 설치 (RTX 3090 / cuEquivariance 호환)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# ✅ MACE 본체 설치
pip install mace-torch nglview ipywidgets rdkit x3dase
pip install -U numpy==2.0
# ✅ cuEquivariance GPU 가속 설치
pip install cuequivariance cuequivariance-torch cuequivariance-ops-torch-cu12

# ✅ 기타 도구
pip install nglview ipywidgets x3dase
conda install -c conda-forge rdkit -y


python -c "import mace; print('✅ MACE installed OK!')"

