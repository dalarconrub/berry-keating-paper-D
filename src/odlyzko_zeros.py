"""
Script para descargar y leer los datasets de zeros de Odlyzko.

Datasets disponibles:
  zeros1     : primeros 100,000 zeros (N=1..100k),      precisión ~3e-9
  zeros6     : primeros 2,001,052 zeros (N=1..2M),      precisión ~4e-9
  zeros3     : 10,000 zeros cerca de N=10^12  → logT ≈ 26.1  (T ≈ 2.67e11)
  zeros4     : 10,000 zeros cerca de N=10^21  → logT ≈ 48.4
  zeros5     : 10,000 zeros cerca de N=10^22  → logT ≈ 51.0

Limitación: zeros3/4/5 solo tienen ~10,000 zeros por bloque.
Si necesitas 200k zeros a esas alturas, no existe fuente pública conocida.
"""

import gzip
import math
import urllib.request
from pathlib import Path

BASE_URL  = "https://www-users.cse.umn.edu/~odlyzko/zeta_tables/"
SAVE_DIR  = Path("odlyzko_data")
SAVE_DIR.mkdir(exist_ok=True)

# Metadatos de cada fichero
# base_N  : índice del primer zero en el fichero
# base_T  : valor que se suma a cada línea para obtener γ (0 si el fichero ya tiene γ completo)
# n_zeros : número aproximado de zeros en el fichero
DATASETS = {
    "zeros1": {
        "url":      BASE_URL + "zeros1",
        "gz_url":   BASE_URL + "zeros1.gz",
        "base_N":   1,
        "base_T":   0.0,          # el fichero contiene γ directamente
        "n_zeros":  100_000,
        "desc":     "Primeros 100,000 zeros",
    },
    "zeros6": {
        "url":      BASE_URL + "zeros6",
        "gz_url":   BASE_URL + "zeros6.gz",
        "base_N":   1,
        "base_T":   0.0,
        "n_zeros":  2_001_052,
        "desc":     "Primeros 2,001,052 zeros",
    },
    "zeros3": {
        "url":      BASE_URL + "zeros3",
        "gz_url":   None,
        "base_N":   10**12 + 1,
        "base_T":   267_653_395_647.0,   # offset publicado por Odlyzko
        "n_zeros":  10_000,
        "desc":     "10,000 zeros cerca de N=10^12  (logT≈26.1)",
    },
    "zeros4": {
        "url":      BASE_URL + "zeros4",
        "gz_url":   None,
        "base_N":   10**21 + 1,
        "base_T":   1_370_919_909_931_995_308_226.8,
        "n_zeros":  10_000,
        "desc":     "10,000 zeros cerca de N=10^21 (logT≈48.4)",
    },
    "zeros5": {
        "url":      BASE_URL + "zeros5",
        "gz_url":   None,
        "base_N":   10**22 + 1,
        "base_T":   15_202_440_115_920_418_066_994.0,
        "n_zeros":  10_000,
        "desc":     "10,000 zeros cerca de N=10^22 (logT≈51.0)",
    },
}


# ── Cobertura logT ─────────────────────────────────────────────────────────
def _T_from_N(N):
    """Aproximación iterativa de T dado N (fórmula de Backlund)."""
    if N < 10:
        return 14.0
    T = 2 * math.pi * N / math.log(N)
    for _ in range(15):
        T = 2 * math.pi * N / math.log(T / (2 * math.pi * math.e))
    return T


def info():
    """Muestra cobertura de cada dataset."""
    print(f"{'Dataset':<10} {'N_inicio':>15} {'T_aprox':>14} {'logT':>6}  Descripción")
    print("─" * 80)
    for name, d in DATASETS.items():
        N  = d["base_N"]
        T  = d["base_T"] if d["base_T"] > 0 else _T_from_N(N)
        lT = math.log(T) if T > 0 else 0
        print(f"  {name:<8} {N:>15,} {T:>14.3e} {lT:>6.2f}  {d['desc']}")


# ── Descarga ───────────────────────────────────────────────────────────────
def _download(url, dest: Path):
    print(f"  Descargando {dest.name} ...")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as r, open(dest, "wb") as f:
            total = 0
            while chunk := r.read(65536):
                f.write(chunk)
                total += len(chunk)
                if total % (10 * 1024 * 1024) == 0:
                    print(f"    {total / 1024 / 1024:.0f} MB...")
        print(f"  ✓ {dest.name}  ({dest.stat().st_size / 1024 / 1024:.1f} MB)")
        return True
    except Exception as e:
        print(f"  ✗ Error: {e}")
        if dest.exists():
            dest.unlink()
        return False


def download(names=None, prefer_gz=True):
    """
    Descarga los datasets indicados (o todos si names=None).
    Si prefer_gz=True intenta la versión comprimida primero.
    """
    if names is None:
        names = list(DATASETS)
    for name in names:
        d    = DATASETS[name]
        dest = SAVE_DIR / name
        gz   = SAVE_DIR / (name + ".gz")

        if dest.exists():
            print(f"  Ya existe: {dest.name}  ({dest.stat().st_size / 1024 / 1024:.1f} MB)")
            continue

        ok = False
        if prefer_gz and d["gz_url"]:
            ok = _download(d["gz_url"], gz)
            if ok:
                print(f"  Descomprimiendo {gz.name} ...")
                with gzip.open(gz, "rb") as fin, open(dest, "wb") as fout:
                    fout.write(fin.read())
                gz.unlink()
                print(f"  ✓ {dest.name}  ({dest.stat().st_size / 1024 / 1024:.1f} MB)")
                ok = True

        if not ok:
            _download(d["url"], dest)


# ── Lectura ────────────────────────────────────────────────────────────────
def read_zeros(name, max_zeros=None):
    """
    Genera pares (N, γ) para el dataset 'name'.
    γ es el valor real de la parte imaginaria del zero (base_T + valor_leído).
    """
    d    = DATASETS[name]
    path = SAVE_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"{path} no encontrado — ejecuta download(['{name}']) primero")

    base_T = d["base_T"]
    N      = d["base_N"]
    count  = 0

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            gamma = base_T + float(line)
            yield N, gamma
            N     += 1
            count += 1
            if max_zeros and count >= max_zeros:
                break


# ── CLI ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys

    cmd = sys.argv[1] if len(sys.argv) > 1 else "info"

    if cmd == "info":
        info()

    elif cmd == "download":
        names = sys.argv[2:] if len(sys.argv) > 2 else None
        download(names)

    elif cmd == "read":
        name = sys.argv[2] if len(sys.argv) > 2 else "zeros3"
        n    = int(sys.argv[3]) if len(sys.argv) > 3 else 10
        for N, gamma in read_zeros(name, max_zeros=n):
            logT = math.log(gamma)
            print(f"  N={N:>25,}  γ={gamma:.10f}  logT={logT:.4f}")

    else:
        print(f"Uso: {sys.argv[0]} [info | download [names...] | read <name> [n]]")
