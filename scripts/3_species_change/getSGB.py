import pickle
import bz2
import subprocess
import shutil
import glob

DB_BASE = "db/mpa_vJun23_CHOCOPhlAnSGB_202403"
SGB = "SGB9340"
OUT_FASTA = f"{SGB}_markers.fna"

# Find bowtie2-inspect: check PATH first, then common conda/pixi locations
bt2_inspect = shutil.which("bowtie2-inspect")
if bt2_inspect is None:
    candidates = glob.glob(
        "/Users/barbaraverhaar/*/envs/*/bin/bowtie2-inspect"
    ) + glob.glob("/opt/conda/envs/*/bin/bowtie2-inspect")
    if candidates:
        bt2_inspect = candidates[0]

if bt2_inspect is None:
    raise RuntimeError(
        "bowtie2-inspect not found. Install bowtie2 or activate the metaphlan environment."
    )

print(f"Using: {bt2_inspect}")

# Get marker IDs for this SGB from the pkl
with bz2.open(f"{DB_BASE}.pkl", "rb") as f:
    db = pickle.load(f)

sgb_markers = [marker for marker in db['markers'] if SGB in marker]
print(f"Found {len(sgb_markers)} markers for {SGB}")

# Stream full index and filter for SGB48820 markers
sgb_markers_set = set(sgb_markers)
written = 0

proc = subprocess.Popen(
    [bt2_inspect, DB_BASE],
    stdout=subprocess.PIPE, text=True
)

with open(OUT_FASTA, "w") as out:
    write = False
    for line in proc.stdout:
        if line.startswith(">"):
            write = SGB in line
            if write:
                written += 1
        if write:
            out.write(line)

proc.wait()
print(f"Written {written} sequences to {OUT_FASTA}")
