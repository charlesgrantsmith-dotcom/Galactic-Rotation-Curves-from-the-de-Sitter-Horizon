# Galactic Rotation Curves from the de Sitter Horizon Temperature

10.5281/zenodo.18806314 

The MOND acceleration constant a₀ ≈ 1.2 × 10⁻¹⁰ m/s² is numerically close to cH₀/(2π) — the acceleration associated with the Gibbons–Hawking temperature of the cosmological horizon. This repository tests that coincidence quantitatively against 171 galaxies from the SPARC database and 14 high-redshift galaxies from Genzel et al. (2017, 2020) and Übler et al. (2017).

## Key Result

The derived constant a_u = cH₀/(2π) = 1.04 × 10⁻¹⁰ m/s² **outperforms** the empirical a₀ = 1.2 × 10⁻¹⁰ m/s² — including within MOND's own interpolation function.

### Fit quality (171 SPARC galaxies)

| Configuration | Med χ²/ν | Global free params |
|---|---|---|
| **MOND simple + a_u** | **1.41** | **0** |
| **Holographic + a_u** | **1.45** | **0** |
| MOND standard + a₀ | 1.52 | 1 |
| Holographic + a₀ | 1.53 | 1 |
| MOND simple + a₀ | 1.53 | 1 |

### 2×2 Isolation test

| Comparison | Winner | Win rate |
|---|---|---|
| a_u vs a₀ in MOND simple (same function) | a_u | 102/171 (59.6%) |
| a_u vs a₀ in Holographic (same function) | a_u | 104/171 (60.8%) |
| Holographic vs MOND at a₀ (same scale) | Holographic | 113/171 (66.1%) |
| Holographic vs MOND at a_u (same scale) | Holographic | 109/171 (63.7%) |

Both the acceleration scale and the interpolation function independently improve fits. The scale improvement is demonstrable within MOND's own functional form.

### Robustness

The result holds across all error treatments, including raw observational errors with no floor applied (108/171, 63.2% win rate).

### High-redshift test

14 galaxies at z = 0.9–2.5 from Genzel et al. and Übler et al. show declining outer rotation curves. The a_u(z) framework predicts this: stronger acceleration at high z shrinks the transition radius, leaving the outer galaxy in the Keplerian decline regime. Observed turnover radii cluster tightly around the predicted transition radius r_t(z) = √(GM_b/a_u(z)), with median r_turn/r_t = 0.80.

Standard MOND (constant a₀) predicts the opposite — flat rotation at *larger* radii at high z.


## Repository structure

```
sparc-gibbons-hawking/
├── README.md
├── LICENSE
├── paper/
│   └── manuscript.md             # Full paper (markdown)
├── scripts/
│   ├── sparc_holographic_fit.py  # Primary 171-galaxy analysis
│   ├── isolation_test.py         # 2×2 scale vs function separation
│   ├── robustness.py             # 5 error treatment tests
│   └── high_z_analysis.py        # High-z rotation curve test
└── results/
    ├── fit_summary.csv           # Per-galaxy results (all models)
    ├── global_statistics.txt     # Summary statistics
    ├── highz_summary.txt         # High-z analysis summary
    ├── sparc_figures/            # 7 SPARC figures (PNG + PDF)
    └── highz_figures/            # 5 high-z figures (PNG + PDF)
```


## Requirements

```bash
pip install numpy scipy matplotlib
```


## Data

Download SPARC mass models from http://astroweb.cwru.edu/SPARC/

Required file: `MassModels_Lelli2016c.mrt` (Newtonian Mass Models, Table 2)

Place in a `SPARC/` directory alongside the scripts, or edit the file path in each script.


## Usage

### Primary analysis (171 galaxies, 3 models)
```bash
python scripts/sparc_holographic_fit.py
```

### 2×2 isolation test (scale vs. function)
```bash
python scripts/isolation_test.py
```

### Robustness check (5 error treatments)
```bash
python scripts/robustness.py
```

### High-redshift analysis
```bash
python scripts/high_z_analysis.py
```


## The physics

In de Sitter spacetime, the cosmological horizon has a Gibbons–Hawking temperature T_GH = ℏH/(2πk_B). The associated acceleration is:

**a_u = k_B T_GH / (ℏ/c) = cH₀/(2π) ≈ 1.04 × 10⁻¹⁰ m/s²**

This is not fitted. It is derived from H₀ = 67.4 km/s/Mpc alone.

The holographic interpolation function:

**g_eff = g_N / (1 − exp(−√(|g_N|/a_u)))**

has the same Newtonian and deep-MOND limits as standard MOND but transitions more sharply. The 2×2 test shows that both the scale (a_u vs a₀) and the function (holographic vs MOND simple) independently improve fits.


## Falsifiable prediction

**a_u(z) = cH(z)/(2π)**

At z = 1: a_u ≈ 1.8× local. At z = 2: a_u ≈ 3.0× local.

This predicts:
- Transition radius shrinks at high z: r_t(z) = √(GM_b/a_u(z))
- Outer rotation curves decline at high z (observed by Genzel et al.)
- BTFR zero-point evolves with H(z)

MOND predicts no evolution (a₀ = constant). CDM predicts halo-dependent evolution unrelated to H(z).


## Citation

If you use this code or results, please cite the SPARC database:
- Lelli, McGaugh & Schombert, Astron. J. 152, 157 (2016)


## License

MIT
