#!/usr/bin/env python3
"""
E4 c parcial: integra c = ∫∫ r × p₂^GUE × (b/R₃) ds₁ds₂
usando los puntos ancla multi-L como fuente de b(v₁,v₂).

Ejecutar:
  python scripts/e4_c_parcial_desde_anclas.py

Lee: scripts/e4_anclas_multi_L_results.json
"""
import sys, os, json
import numpy as np
from numpy.polynomial.legendre import leggauss

# ── p₂^GUE exacta (Schur complement) ────────────────────────────────────────
def sine_kernel(x, y):
    r = x - y
    return 1.0 if abs(r) < 1e-12 else np.sin(np.pi * r) / (np.pi * r)

def bornemann_det(kernel, a, b_lim, n_quad=32):
    xi, wi = leggauss(n_quad)
    mid, half = (a + b_lim) / 2, (b_lim - a) / 2
    t, w = mid + half * xi, half * wi
    K = np.array([[kernel(t[i], t[j]) * np.sqrt(w[i] * w[j])
                   for j in range(n_quad)] for i in range(n_quad)])
    return np.linalg.det(np.eye(n_quad) - K)

def p2_GUE(s1, s2, n_quad=32):
    if s1 < 1e-6 or s2 < 1e-6:
        return 0.0
    S = s1 + s2
    pts = [0.0, s1, S]
    M3 = np.array([[sine_kernel(pts[i], pts[j]) for j in range(3)] for i in range(3)])
    det3 = np.linalg.det(M3)
    if det3 < 0:
        return 0.0
    try:
        M3_inv = np.linalg.inv(M3)
    except np.linalg.LinAlgError:
        return 0.0
    def K_cond(x, y):
        kx = np.array([sine_kernel(x, p) for p in pts])
        ky = np.array([sine_kernel(p, y) for p in pts])
        return sine_kernel(x, y) - kx @ M3_inv @ ky
    return det3 * bornemann_det(K_cond, 0.0, S, n_quad)

def R3_GUE(v1, v2):
    S = np.sinc
    s01 = S(v1); s02 = S(v2); s12 = S(v1 - v2)
    return 1 - s01**2 - s02**2 - s12**2 + 2 * s01 * s02 * s12

# ── Cargar anclas ────────────────────────────────────────────────────────────
ANCHOR_FILE = os.path.join(os.path.dirname(__file__), 'e4_anclas_multi_L_results.json')

with open(ANCHOR_FILE) as f:
    anchors = json.load(f)

print(f'Anclas cargados: {len(anchors)} puntos')
print()

# ── Interpolador b/R₃ por distancia inversa ──────────────────────────────────
pts_coords = np.array([[r['v1'], r['v2']] for r in anchors])
pts_bR3 = np.array([r['fit_3p']['b'] / r['R3_GUE'] for r in anchors])

def b_over_R3_interp(v1, v2):
    dists = np.sqrt((pts_coords[:, 0] - v1)**2 + (pts_coords[:, 1] - v2)**2)
    w = 1.0 / (dists + 0.3)**2
    return np.sum(w * pts_bR3) / np.sum(w)

# ── Tabla de anclas ──────────────────────────────────────────────────────────
print('%-12s %8s %10s %10s' % ('(v1,v2)', 'R3_GUE', 'b_3p', 'b/R3'))
print('-' * 44)
for r in sorted(anchors, key=lambda x: (x['v1'], x['v2'])):
    f3 = r['fit_3p']
    print('(%.1f,%.1f)   %8.4f %+10.2f %+10.2f' %
          (r['v1'], r['v2'], r['R3_GUE'], f3['b'], f3['b'] / r['R3_GUE']))
print()

# ── Integral 2D ──────────────────────────────────────────────────────────────
N_INT = 20
s_max = 3.5
R3_THRESH = 0.15
c_emp = 1.2446

ds = s_max / N_INT
s_grid = np.linspace(ds / 2, s_max - ds / 2, N_INT)

print(f'Integrando c = integral r * p2 * (b/R3) ds1 ds2')
print(f'N_INT={N_INT}, s_max={s_max}, R3_thresh={R3_THRESH}')
print()

c_sum = 0.0
norm_sum = 0.0
n_used = 0

for i, s1 in enumerate(s_grid):
    for j, s2 in enumerate(s_grid):
        if s1 < 0.1 or s2 < 0.1:
            continue
        v1, v2 = s1, s1 + s2
        R3g = R3_GUE(v1, v2)
        if R3g < R3_THRESH:
            continue
        r_val = min(s1, s2) / max(s1, s2)
        p2_val = p2_GUE(s1, s2, n_quad=24)
        bR3_val = b_over_R3_interp(v1, v2)
        c_sum += r_val * p2_val * bR3_val * ds**2
        norm_sum += r_val * p2_val * ds**2
        n_used += 1

print('=' * 55)
print(f'RESULTADO: c parcial ({len(anchors)} anclas multi-L)')
print('=' * 55)
print(f'  c_parcial   = {c_sum:+.4f}')
print(f'  c_emp       =  {c_emp:.4f} +- 0.040')
print(f'  ratio       =  {c_sum / c_emp:.3f}  ({c_sum / c_emp * 100:.1f}% de c_emp)')
print(f'  norm (r*p2) =  {norm_sum:.4f}')
print(f'  puntos      =  {n_used}/{N_INT**2}')
print('=' * 55)
print()

# ── Contribucion por ancla ───────────────────────────────────────────────────
print('Contribucion por ancla (peso en la integral):')
for r in sorted(anchors, key=lambda x: (x['v1'], x['v2'])):
    v1a, v2a = r['v1'], r['v2']
    bR3 = r['fit_3p']['b'] / r['R3_GUE']
    s1a = v1a
    s2a = v2a - v1a
    if s2a > 0.1:
        p2a = p2_GUE(s1a, s2a, n_quad=24)
        ra = min(s1a, s2a) / max(s1a, s2a)
        contrib = ra * p2a * bR3
        print(f'  ({v1a:.1f},{v2a:.1f}): b/R3={bR3:+8.2f}  p2={p2a:.3f}  '
              f'r={ra:.3f}  r*p2*b/R3={contrib:+.3f}')

# ── Estimacion simple ────────────────────────────────────────────────────────
bR3_mean = np.mean(pts_bR3)
print()
print(f'Estimacion simple: c ~ mean(b/R3) * integral(r*p2)')
print(f'  = {bR3_mean:.2f} * {norm_sum:.4f} = {bR3_mean * norm_sum:.4f}')
