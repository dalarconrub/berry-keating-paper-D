"""
02 — Mechanism of the correction (Paper D, Section 3)

Decomposition c = c_std + c_corr = +1.60 + (-0.36) = 1.24.
Same as Paper A Section 5 but summarized.

Run from repo root: python src/02_mechanism.py
"""

dr_dstd = -0.749;  a_std = -2.13
dr_dcorr = +0.116; a_corr = -3.08

c_std = dr_dstd * a_std
c_corr = dr_dcorr * a_corr
c_total = c_std + c_corr

print("=== Decomposition of c ===")
print(f"c_std  = ({dr_dstd}) x ({a_std}) = {c_std:+.3f} (+128%)")
print(f"c_corr = ({dr_dcorr}) x ({a_corr}) = {c_corr:+.3f} (-29%)")
print(f"c_total = {c_total:.3f} (99.5% of c_emp = 1.245)")
