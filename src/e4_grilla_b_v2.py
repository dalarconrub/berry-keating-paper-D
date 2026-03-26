#!/usr/bin/env python3
"""
E4 GRILLA v2: b(v₁,v₂) con mejor condicionamiento para c ab initio.

Mejoras respecto a v1:
  - L_vals = [18, 22, 26, 30, 34]   (ratio 34/18 = 1.89; todos L altos → menos ruido)
  - v_grid denso: paso 0.1 en [0.3, 2.5]  (23 puntos vs 10 → ~450 puntos válidos)
  - Richardson implícito: pares (18,26), (22,30) con ratio ≈ √2 para eliminar a/L

Modelo ajustado: δR₃ = a/L + b/L²
Salida: scripts/e4_grilla_b_v2_results.json

TIEMPO ESTIMADO: ~6-8 horas. Ejecutar overnight:
  python scripts/e4_grilla_b_v2.py
"""
import mpmath, numpy as np, json, time, sys
from scipy.optimize import curve_fit

mpmath.mp.dps = 20

def sieve(N):
    is_p = np.ones(N+1, dtype=bool); is_p[0]=is_p[1]=False
    for i in range(2, int(N**0.5)+1):
        if is_p[i]: is_p[i*i::i] = False
    return np.where(is_p)[0]

PRIMES = sieve(50000).astype(float)
NP = 3000

def A_f(x):
    r = mpmath.mpf(1)
    for p in PRIMES[:NP]:
        p = mpmath.mpf(p); p1x = mpmath.power(p, 1+x)
        r *= (1 - 1/p1x) * (1 - 2/p + 1/p1x) / (1 - 1/p)**2
    return r

def B_f(x):
    r = mpmath.mpf(0)
    for p in PRIMES[:NP]:
        p = mpmath.mpf(p)
        r += (mpmath.log(p) / (mpmath.power(p, 1+x) - 1))**2
    return r

def Q_f(x, y):
    r = mpmath.mpf(0)
    for p in PRIMES[:NP]:
        p = mpmath.mpf(p)
        r -= mpmath.log(p)**3 / (
            mpmath.power(p, 2+x+y) *
            (1 - 1/mpmath.power(p, 1+x)) *
            (1 - 1/mpmath.power(p, 1+y))
        )
    return r

def P_f(x, y):
    Ax = A_f(x); S = mpmath.mpf(0)
    for p in PRIMES[:NP]:
        p = mpmath.mpf(p); lp = mpmath.log(p)
        px = mpmath.power(p, x); py = mpmath.power(p, y)
        p1x = mpmath.power(p, 1+x); p1y = mpmath.power(p, 1+y)
        num = (1 - 1/px) * (1 - 1/px - 1/py + 1/p1y) * (-lp)
        den = (1 - 1/mpmath.power(p, 1-x+y)) * (1 - 1/p1y) * \
              (1 - 2/p + 1/p1x) * mpmath.power(p, 2-x+y)
        if abs(den) > 1e-30:
            S += num / den
    return Ax * S

def zpz(s):
    return mpmath.zeta(s, derivative=1) / mpmath.zeta(s)

def zpz_prime(s):
    z = mpmath.zeta(s)
    zp = mpmath.zeta(s, derivative=1)
    zpp = mpmath.zeta(s, derivative=2)
    return zpp/z - (zp/z)**2

def I_closed(a1, a2, beta, T):
    L = mpmath.log(T/(2*mpmath.pi)); s1 = a1+beta; s2 = a2+beta; r = mpmath.mpc(0)
    for sa, sb, aa, ab in [(s1,s2,a1,a2), (s2,s1,a2,a1)]:
        if abs(sa) < 1e-20: continue
        ph = T * mpmath.exp(-sa*L) / (1-sa)
        zz = mpmath.zeta(1-sa) * mpmath.zeta(1+sa)
        Qv = Q_f(sa, sb); Av = A_f(sa); Pv = P_f(sa, sb)
        if abs(ab-aa) > 1e-20 and abs(ab+beta) > 1e-20:
            zd = zpz(1+ab-aa) - zpz(1+ab+beta)
        elif abs(ab-aa) > 1e-20:
            zd = zpz(1+ab-aa) + 1/(ab+beta)
        else:
            zd = mpmath.mpf(0)
        r += ph * zz * (Qv + Av*zd + Pv)
    return r

def I1_corrected(alpha, beta, T):
    L = mpmath.log(T/(2*mpmath.pi)); s = alpha+beta
    if abs(s) < 1e-20: s = mpmath.mpf('1e-10')
    zpzp = zpz_prime(1+s)
    int_log = T * (L-1)
    exp_sL = mpmath.exp(-s*L)
    int_log_phase = T * exp_sL * (L/(1-s) + 1/(1-s)**2)
    zz = mpmath.zeta(1+s) * mpmath.zeta(1-s)
    bracket2 = zz * A_f(s) - B_f(s)
    return int_log * zpzp + int_log_phase * bracket2

def R3_GUE(v1, v2):
    S = np.sinc
    return 1 - S(v1)**2 - S(v2)**2 - S(v1-v2)**2 + 2*S(v1)*S(v2)*S(v1-v2)

def R3_CS(v1, v2, Lv):
    """R₃ de Conrey-Snaith para (v₁,v₂) y L dado."""
    Tv = mpmath.exp(Lv) * 2 * mpmath.pi
    rho = Lv / (2*np.pi)
    a1 = mpmath.mpc(0, 2*mpmath.pi*v1/Lv)
    a2 = mpmath.mpc(0, 2*mpmath.pi*v2/Lv)
    z  = mpmath.mpf(0)
    Is = mpmath.re(
        I_closed(a1, a2, z, Tv) + I_closed(z, a1, -a2, Tv) + I_closed(z, a2, -a1, Tv) +
        I_closed(-a1, -a2, z, Tv) + I_closed(z, -a2, a1, Tv) + I_closed(z, -a1, a2, Tv)
    )
    I1s = mpmath.re(
        I1_corrected(z, a2, Tv) + I1_corrected(z, a1, Tv) + I1_corrected(-a2, a1, Tv) +
        I1_corrected(-a2, z, Tv) + I1_corrected(-a1, a2, Tv) + I1_corrected(-a1, z, Tv)
    )
    log3 = float(Tv) * (Lv**3 - 3*Lv**2 + 6*Lv - 6)
    cont = log3 + float(Is) + float(I1s)
    return cont / ((2*np.pi)**3 * float(Tv) * rho**3)

# ============================================================
# Configuración de grilla v2
# ============================================================
v_grid = list(np.round(np.arange(0.3, 2.6, 0.1), 2))   # 23 puntos: 0.3, 0.4, ..., 2.5
L_vals = [18, 22, 26, 30, 34]                             # ratio 34/18 = 1.89
OUTPUT = 'scripts/e4_grilla_b_v2_results.json'

# Filtrar pares válidos
pairs = []
for v1 in v_grid:
    for v2 in v_grid:
        if abs(v1-v2) < 0.2 or v1 < 0.2 or v2 < 0.2:
            continue
        pairs.append((v1, v2))

N_pairs = len(pairs)
N_L = len(L_vals)
t_est = N_pairs * N_L * 3 / 60   # ~3s por evaluación

print(f'Grilla v2: {len(v_grid)} valores de v, paso 0.1')
print(f'Pares válidos: {N_pairs}')
print(f'L_vals = {L_vals}  (ratio {L_vals[-1]/L_vals[0]:.2f})')
print(f'Estimado: {t_est:.0f} minutos (~{t_est/60:.1f} horas)')
print('='*60)
sys.stdout.flush()

# Cargar progreso previo si existe
try:
    with open(OUTPUT, encoding='utf-8') as f:
        saved = json.load(f)
    results = saved.get('results', {})
    print(f'Retomando: {len(results)} puntos ya calculados')
except FileNotFoundError:
    results = {}

t_total = time.time()
done = 0

for idx, (v1, v2) in enumerate(pairs):
    key = f'{v1:.1f},{v2:.1f}'
    if key in results:
        done += 1
        continue

    t0 = time.time()
    R3g = R3_GUE(v1, v2)

    R3_vals = []
    for Lv in L_vals:
        R3 = R3_CS(v1, v2, Lv)
        R3_vals.append(R3)

    # Ajuste δR₃ = a/L + b/L²
    Ls = np.array(L_vals, dtype=float)
    dR3 = np.array(R3_vals) - R3g

    try:
        def m12(L, a, b): return a/L + b/L**2
        p, cov = curve_fit(m12, Ls, dR3)
        a_coeff, b_coeff = p[0], p[1]
        a_err, b_err = np.sqrt(np.diag(cov))
    except Exception:
        a_coeff, b_coeff = 0.0, 0.0
        a_err, b_err = 99.0, 99.0

    # Richardson: eliminar a/L usando par (L, L·√2)
    # Par (18,26): ratio 26/18 = 1.444 ≈ √2
    # b_rich = [R3(L2) - (L1/L2)·R3(L1) - R3g·(1-L1/L2)] / (1/L1² - 1/L2²·...)
    # Nota: guardamos R3_vals para que el notebook pueda hacer Richardson
    b_rich = None
    # Par índice 0 y 2: L=18 y L=26, ratio=26/18=1.444
    if len(L_vals) >= 3:
        L1, L2 = L_vals[0], L_vals[2]   # 18, 26
        r1, r2 = R3_vals[0], R3_vals[2]
        # Eliminar término a/L: b = [dR3(L2)·L2² - dR3(L1)·L1²·(L2/L1)] / (L2² - L1²·(L2/L1))
        # Simplificando con dR3=a/L+b/L²:
        # dR3(L1)·L1 - dR3(L2)·L2 = b·(1/L1 - 1/L2) → b = [dR3(L1)·L1 - dR3(L2)·L2]/(1/L1-1/L2)
        dr1 = r1 - R3g; dr2 = r2 - R3g
        denom = 1/L1 - 1/L2
        if abs(denom) > 1e-10:
            b_rich = (dr1*L1 - dr2*L2) / (L1*L2*denom)
            # Equivalente: b_rich = (dr1*L1² - dr2*L2²) · L1·L2 / (L2²-L1²) ... simplificado:
            b_rich = (dr1*L1 - dr2*L2) * L1*L2 / (L2 - L1)

    dt = time.time() - t0
    done += 1

    results[key] = {
        'v1': v1, 'v2': v2, 'R3_GUE': R3g,
        'a': a_coeff, 'b': b_coeff,
        'a_err': a_err, 'b_err': b_err,
        'b_rich': b_rich,
        'R3_vals': R3_vals
    }

    elapsed = time.time() - t_total
    remaining = elapsed / done * (N_pairs - done) / 60
    print(f'({v1:.1f},{v2:.1f}): a={a_coeff:+.4f} b={b_coeff:+.3f} '
          f'b_rich={b_rich:+.3f} R3g={R3g:.4f} ({dt:.0f}s) '
          f'[{done}/{N_pairs}, ~{remaining:.0f}min]')
    sys.stdout.flush()

    with open(OUTPUT, 'w', encoding='utf-8') as f:
        json.dump({'v_grid': v_grid, 'L_vals': L_vals, 'results': results}, f, indent=2)

dt_total = time.time() - t_total
print(f'\nTotal: {dt_total/60:.0f} min')
print(f'Guardado en {OUTPUT}')

# ── Estadísticas finales ──────────────────────────────────────────────────────
print(f'\n{"="*60}')
print('Estadísticas b(v₁,v₂):')
bs = [r['b'] for r in results.values()]
bs_rich = [r['b_rich'] for r in results.values() if r['b_rich'] is not None]
print(f'  N puntos:       {len(bs)}')
print(f'  mean(b):        {np.mean(bs):+.3f}')
print(f'  std(b):         {np.std(bs):.3f}')
print(f'  mean(b_rich):   {np.mean(bs_rich):+.3f}')
print(f'  std(b_rich):    {np.std(bs_rich):.3f}')

# Región física clave: v₁∈[0.7,1.3], v₂∈[1.4,2.6] (gaps típicos GUE)
bs_phys = [r['b'] for r in results.values()
           if 0.7 <= r['v1'] <= 1.3 and 1.4 <= r['v2'] <= 2.6]
print(f'  mean(b) región física [0.7–1.3, 1.4–2.6]: {np.mean(bs_phys):+.3f}  ({len(bs_phys)} pts)')
print('='*60)
