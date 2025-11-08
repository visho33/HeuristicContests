import os
import random
import subprocess
from pathlib import Path
from subprocess import TimeoutExpired

# Configuración
TRAIN_FOLDERS = ["training-01", "training-02"]
OUTPUT_IN_DIR = Path("train/in")
OUTPUT_GT_DIR = Path("train/gt")
GENERATE_SCRIPT = "generate_data.py"
GT_SCRIPT = "ground_truth_generator.py"
TIMEOUT = 600  # tiempo máximo (segundos) por bloque (para .in y .gt)

# Aseguramos que las carpetas existan
OUTPUT_IN_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_GT_DIR.mkdir(parents=True, exist_ok=True)

def get_max_offset_seconds(file_path):
    """Obtiene el offset máximo (en segundos) de un feed."""
    with open(file_path, "rb") as f:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b"\n":
            f.seek(-2, os.SEEK_CUR)
        last_line = f.readline().decode().strip()
    offset_micro = int(last_line.split()[0])
    return offset_micro / 1_000_000

def pick_random_feed_and_skip():
    """Elige un feed al azar y un skip_s válido según el offset máximo."""
    folder = random.choice(TRAIN_FOLDERS)
    candidates = [
        str(Path(folder) / f)
        for f in os.listdir(folder)
        if f.endswith(".txt")
    ]
    feed_file = random.choice(candidates)
    max_offset_sec = get_max_offset_seconds(feed_file)

    # Margen: 900 (bloque) + 1600 (para horizons y ground truth)
    max_skip = int(max_offset_sec) - (900 + 1600)
    if max_skip <= 0:
        raise ValueError(f"Archivo {feed_file} demasiado corto para un bloque válido.")
    skip_s = random.randint(0, max_skip)
    return feed_file, skip_s

def try_generate_pair(idx, feed_file, skip_s):
    """Intenta generar un par .in/.gt, devuelve True si se completó correctamente."""
    out_in = OUTPUT_IN_DIR / f"{idx}.in"
    out_gt = OUTPUT_GT_DIR / f"{idx}.gt"

    # Eliminar si quedaron de intentos previos
    if out_in.exists():
        out_in.unlink()
    if out_gt.exists():
        out_gt.unlink()

    try:
        # Generar .in con timeout
        subprocess.run([
            "python3", GENERATE_SCRIPT,
            "--input", feed_file,
            "--output", str(out_in),
            "--skip_s", str(skip_s)
        ], check=True, timeout=TIMEOUT)

        # Generar .gt con timeout
        subprocess.run([
            "python3", GT_SCRIPT,
            str(feed_file),
            str(out_gt),
            "--skip_s", str(skip_s)
        ], check=True, timeout=TIMEOUT)

        return True  # todo salió bien

    except (TimeoutExpired, subprocess.CalledProcessError) as e:
        print(f"  ⚠️  Falló bloque {idx} (feed={feed_file}, skip={skip_s}): {e}")
        # Limpiar archivos incompletos
        if out_in.exists():
            out_in.unlink()
        if out_gt.exists():
            out_gt.unlink()
        return False

def make_train_dataset(n_files):
    """Genera exactamente n_files pares .in/.gt (reintenta hasta lograrlo)."""
    i = 265
    attempts = 0
    while i < n_files:
        attempts += 1
        feed_file, skip_s = pick_random_feed_and_skip()
        idx = f"{i:04d}"
        print(f"[{i+1}/{n_files}] Intento #{attempts} | Feed: {feed_file} | skip_s={skip_s}")

        if try_generate_pair(idx, feed_file, skip_s):
            i += 1  # sólo avanzamos si se completó correctamente
        else:
            print(f"  → Reintentando con otro feed y skip (no se cuenta para el total).")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Genera un dataset train/ con archivos .in y .gt")
    parser.add_argument("--n", type=int, required=True, help="Número de archivos a generar")
    args = parser.parse_args()

    random.seed(33)
    make_train_dataset(args.n)
