# Berry-Keating convergence implies the Riemann Hypothesis

**Paper:** "Empirical proof that Berry-Keating convergence implies the Riemann Hypothesis"

**Author:** David Alarcon, Universidad Pablo de Olavide, Sevilla, Spain

**Submitted to:** Nature

## Abstract

We report the first precision measurement of the rate at which the gap ratio
statistic of Riemann zeta zeros converges to the GUE prediction. Using
high-precision zeros up to height T ~ 3x10^10, we find the convergence rate
c/log^2(T) with c = 1.245 +/- 0.040. We identify the mechanism (narrowing of
the spacing distribution) and prove that if this convergence pattern holds
exactly, the Riemann Hypothesis is true.

## Main results

1. **First measurement:** c = 1.245 +/- 0.040, R_inf 6.1 sigma below R_GUE
2. **Mechanism:** c_std + c_corr = 1.60 - 0.36 = 1.24 (99.5% of c_emp)
3. **Theorem:** (H1) => RH, using VK + Simon + Cauchy (no open conjectures)
4. **Ab initio:** c = 1.23 from Conrey-Snaith R_3 (99% of c_emp)

## Repository structure

```
├── main.tex              Main manuscript
├── methods.tex           Methods section
├── supplementary.tex     Supplementary Information
├── references.bib        Bibliography (20 references)
├── naturemag.bst         Nature bibliography style
├── cover_letter.tex      Cover letter to editors
├── main.pdf              Compiled manuscript (7 pages)
├── figures/              All figures (PDF + PNG)
│   ├── fig1_r_vs_logT           r(T) with Model A fit
│   ├── fig2_c_decomposition     c = c_std + c_corr
│   ├── fig3_delta_p_s           Narrowing of p(s)
│   ├── fig4_argument_D          Proof diagram D1-D4
│   ├── fig5_exclusion_region    Bounds on off-line fraction f
│   ├── ed_fig1_sigma2_delta3    Number variance and spectral rigidity
│   ├── ed_fig2_rh_bound         RH consistency bound
│   ├── ed_fig3_c_scan_R3        Ab initio c vs R3 threshold
│   └── ed_fig4_std_corr         std(s) and Corr convergence
├── data/
│   ├── dataset_v6_21pts.dat             Dataset (21 points)
│   ├── e4_grilla_b_v2_results.json      Conrey-Snaith grid (436 pairs)
│   └── e4_anclas_multi_L_results.json   Multi-L anchors (29 points)
├── notebooks/                           Analysis code (Jupyter)
│   ├── platt_analisis.ipynb             Main analysis: r(T), model fits
│   ├── c_teorico_bk_directo.ipynb       Mechanism: std, Corr, decomposition
│   ├── b_cero_simetria.ipynb            A1[r]=0 by symmetry (3 proofs)
│   ├── sigma2_delta3_bogomolny.ipynb    Sigma2, Delta3 vs BK saturation
│   ├── tension_Rinf_investigacion.ipynb R_inf < R_GUE diagnostic
│   ├── rh_perturbacion_cota.ipynb       RH consistency bound (eps < 1%)
│   ├── rh_paso7_offline_constraint.ipynb Step 7: off-line zero effect
│   ├── rh_paso8_propagacion.ipynb       Step 8: F(sigma_0) != 0
│   ├── rh_paso9_fraccion_acumulada.ipynb Step 9: bounds on f
│   ├── rh_paso10_teorema_incompatibilidad.ipynb  Step 10: theorem
│   ├── rh_paso11_d2_lipschitz_proof.ipynb Step 11: argument D (proof)
│   ├── rh_baseline_gue.ipynb            GUE baseline (Fredholm, PV)
│   └── e4_R3_empirico_vs_CS.ipynb       Empirical R3 vs Conrey-Snaith
└── src/
    ├── platt_zeros.py           Platt zero reader
    ├── odlyzko_zeros.py         Odlyzko zero reader
    ├── e4_grilla_b_v2.py        Conrey-Snaith R3 grid computation
    └── e4_c_parcial_desde_anclas.py  Ab initio c integration
```

## Compile

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Data

The dataset (`data/dataset_v6_21pts.dat`) contains 21 measurements of the
gap ratio statistic at different heights on the critical line. Columns:
logT, r, sigma, source.

Load in Python:
```python
import numpy as np
logT, r, sigma = np.loadtxt('data/dataset_v6_21pts.dat', usecols=(0,1,2), unpack=True)
```

## Companion papers

- **Paper 2** (Annals): Rigorous version of the theorem with complete proofs
- **Paper 3** (Comm. Math. Phys.): Ab initio derivation of c from Conrey-Snaith R_3

## License

This work is submitted for publication. The dataset is freely available for
academic use.
