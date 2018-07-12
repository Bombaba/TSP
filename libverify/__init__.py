import subprocess
from pathlib import Path

OUT_NAME = "checkdat"

dir_path = Path(__path__[0])
c_path = dir_path / "check_validity_prec.c"
out_path = dir_path / OUT_NAME

process = subprocess.run(["gcc", str(c_path), "-lm", "-o", str(out_path)], check=True)
print("check_validity build success")

def check(tsp_path, dat_path):
    process = subprocess.Popen([str(out_path), str(dat_path), str(tsp_path)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process
