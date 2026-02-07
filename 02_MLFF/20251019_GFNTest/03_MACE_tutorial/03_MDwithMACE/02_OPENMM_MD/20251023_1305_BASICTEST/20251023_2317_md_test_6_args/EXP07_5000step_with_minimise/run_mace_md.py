import torch
import subprocess
import sys

# ✅ autograd 완전 차단
torch.set_grad_enabled(False)

# ✅ 명령 실행
cmd = " ".join(sys.argv[1:])
subprocess.run(cmd, shell=True)

