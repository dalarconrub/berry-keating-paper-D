# Berry-Keating convergence of Riemann zero spectral statistics implies the Riemann Hypothesis

**Author:** David Alarcon, Universidad Pablo de Olavide, Sevilla, Spain

**Venue:** Nature

**Zenodo DOI:** 10.5281/zenodo.19268989

## Abstract

We report the first precision measurement of the rate at which the gap ratio
statistic of Riemann zeta zeros converges to the GUE prediction. Using
high-precision zeros up to height T ~ 3x10^10 (log T = 24), we find
<r>(T) = 0.59891(13) + 1.245(40)/log^2(T). We identify the mechanism
(narrowing of the spacing distribution) and prove that if this convergence
pattern holds exactly, the Riemann Hypothesis is true. In a companion paper,
we show that this convergence follows unconditionally from the Rudnick-Sarnak
theorem and a linear programming bound.

## Repository structure

```
main.tex              Main manuscript (8 pages)
methods.tex           Methods section
supplementary.tex     Supplementary Information (6 pages)
references.bib        Bibliography (19 references)
naturemag.bst         Nature bibliography style
cover_letter.tex      Cover letter to editors
figures/
  fig1_r_vs_logT          <r>(T) with Model A fit (21 points)
  fig2_c_decomposition    c = c_std + c_corr decomposition
  fig4_argument_D         Complete proof structure (Paper C + Theorem 1)
  fig5_exclusion_region   Bounds on off-line zero fraction f
  ed_fig1_sigma2_delta3   Number variance and spectral rigidity
  ed_fig2_rh_bound        RH consistency bound (eps < 1%)
  ed_fig3_c_scan_R3       Two-channel CS decomposition
  ed_fig4_std_corr        std(s) and Corr convergence
data/
  dataset_v6_21pts.dat    Dataset (21 points: logT, r, sigma)
src/
  01_measurement.py       Section 2: dataset and Model A fit
  02_mechanism.py         Section 3: c = c_std + c_corr decomposition
  03_theorem_convergence_implies_rh.py  Section 4: VK + Simon + Cauchy proof
```

## Compile

```bash
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
pdflatex supplementary.tex && pdflatex supplementary.tex
```

## Companion papers

- **Paper A** (Comm. Math. Phys.): Full dataset, model comparison, mechanism analysis
- **Paper B** (J. Number Theory): Ab initio derivation of c from Conrey-Snaith R3
- **Paper C** (PNAS): Unconditional proof: RS correlations + LP imply RH

## License

CC BY 4.0
