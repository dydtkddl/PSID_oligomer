import torch
print("CUDA OK?", torch.cuda.is_available())
import torch
print("Torch version:", torch.__version__)
print("CUDA available:", torch.cuda.is_available())
print("CUDA version torch was built with:", torch.version.cuda)
import cuequivariance_torch
print("cuEquivariance OK")
