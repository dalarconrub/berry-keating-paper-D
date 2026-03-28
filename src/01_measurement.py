"""
01 — Precision measurement (Paper D, Section 2)

Loads dataset v6, fits Model A, displays results.
Same as Paper A Section 2-3 but summarized for the Nature paper.

Run from repo root: python src/01_measurement.py
"""
import numpy as np
from scipy.optimize import curve_fit

data = np.loadtxt("data/dataset_v6_21pts.dat", skiprows=1)
logT, r_mean, sigma = data[:, 0], data[:, 1], data[:, 2]

def model_A(x, Rinf, c):
    return Rinf + c / x**2

popt, pcov = curve_fit(model_A, logT, r_mean, sigma=sigma, absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
chi2 = np.sum(((r_mean - model_A(logT, *popt)) / sigma)**2)
dof = len(logT) - 2

R_GUE = 0.59971
Rinf, c = popt
Rinf_err, c_err = perr

print("=== Model A: <r>(T) = R_inf + c/log^2(T) ===")
print(f"R_inf = {Rinf:.5f} +/- {Rinf_err:.5f}")
print(f"c     = {c:.4f} +/- {c_err:.4f}")
print(f"chi2/dof = {chi2/dof:.3f} ({dof} dof)")
print(f"R_inf - R_GUE = {Rinf - R_GUE:.5f} ({(Rinf - R_GUE)/Rinf_err:.1f} sigma)")
