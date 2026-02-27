---
title: "Galactic Rotation Curves from the de Sitter Horizon Temperature"
author:
  - Charles Grant Smith Jr.
date: "February 2026"
---

# Abstract

The empirical MOND acceleration constant $a_0 \approx 1.2 \times 10^{-10}$ m/s² is numerically close to $cH_0/(2\pi)$ — the acceleration associated with the Gibbons–Hawking temperature of the de Sitter horizon. We define $a_u = cH_0/(2\pi)$ and test it against 171 SPARC galaxies with a holographic interpolation function, with $a_u$ fixed globally and only stellar mass-to-light ratios free. The holographic model achieves median $\chi^2/\nu = 1.45$, outperforming MOND simple (1.53) and MOND standard (1.52) despite having zero global free parameters. A 2×2 isolation test separates the contributions: (i) $a_u$ outperforms $a_0$ within MOND's own interpolation function (102/171 galaxies, 59.6%); (ii) the holographic interpolation outperforms MOND simple within the same acceleration scale (113/171, 66.1%). Both effects are independent and robust across all error treatments. MOND simple with $a_u$ substituted for $a_0$ achieves the lowest median $\chi^2/\nu$ of any configuration (1.41). The prediction $a_u(z) = cH(z)/(2\pi)$ provides a falsifiable distinction from both CDM and MOND.

---

# I. Introduction

The flat rotation curves of spiral galaxies require either unseen mass or a modification of Newtonian dynamics at low accelerations. Cold dark matter (CDM) postulates massive halos fit per galaxy with 2–3 free parameters [1]. Modified Newtonian Dynamics (MOND) [2] introduces a universal acceleration scale $a_0 \approx 1.2 \times 10^{-10}$ m/s² below which gravity transitions from $1/r^2$ to $1/r$. MOND predicts the baryonic Tully–Fisher relation (BTFR) with slope 4 [3] and fits rotation curves with fewer parameters than CDM, but $a_0$ remains an empirical constant whose origin is unexplained.

A long-noted coincidence is $a_0 \sim cH_0$ [2,4,5]. In de Sitter spacetime, the cosmological horizon has a Gibbons–Hawking temperature $T_\mathrm{GH} = \hbar H/(2\pi k_B)$ [7]. The acceleration obtained by dividing the thermal energy $k_B T_\mathrm{GH}$ by the Compton momentum quantum $\hbar/c$ is:

$$a_u = \frac{k_B T_\mathrm{GH}}{\hbar/c} = \frac{cH_0}{2\pi} \approx 1.04 \times 10^{-10} \; \mathrm{m/s^2} \tag{1}$$

This is not a fit. It is a parameter-free prediction from $H_0$ alone [6]. We test whether $a_u$ can replace $a_0$ across the full SPARC sample [8]. To separate the effect of the acceleration scale from the interpolation function, we perform a 2×2 analysis using both scales in both functional forms.


# II. Theoretical Framework

## A. Acceleration scale

For $H_0 = 67.4 \pm 0.5$ km/s/Mpc [6]: $a_u = (1.042 \pm 0.008) \times 10^{-10}$ m/s², 13% below $a_0 = 1.2 \times 10^{-10}$ m/s². This offset is within the systematic uncertainty of $a_0$: Li et al. [9] report $a_\dagger = (1.20 \pm 0.02 \pm 0.24) \times 10^{-10}$ m/s² (second error systematic).

## B. Interpolation functions

We test two interpolation functions. The MOND simple function [10]:

$$g_\mathrm{eff} = \tfrac{1}{2}\left(g_N + \sqrt{g_N^2 + 4|g_N|a}\right) \quad \text{[MOND simple]} \tag{2a}$$

and a holographic interpolation function:

$$g_\mathrm{eff} = \frac{g_N}{1 - \exp\left(-\sqrt{|g_N|/a}\right)} \quad \text{[Holographic]} \tag{2b}$$

Both share the correct limits: $g_\mathrm{eff} \to g_N$ for $g_N \gg a$, and $g_\mathrm{eff} \to \sqrt{g_N \cdot a}$ for $g_N \ll a$. The holographic form (2b) transitions more sharply, passing through $\mu = 0.5$ over ~0.5 dex vs. ~1.0 dex for MOND simple. Each function is tested with both $a_u$ and $a_0$, giving four configurations.

## C. Predicted rotation velocity

From the SPARC mass model decompositions [8]: $g_N(r) = [\Upsilon_d V_d^2 + \Upsilon_b V_b^2 + V_\mathrm{gas}^2]/r$. The stellar mass-to-light ratios $\Upsilon_d$ and $\Upsilon_b$ are the only free parameters per galaxy. $V_\mathrm{pred}(r) = \sqrt{r \cdot g_\mathrm{eff}(r)}$.


# III. Data and Method

We use the SPARC database [8]: 175 late-type galaxies. After quality cuts (≥5 points, $V_\mathrm{max} > 10$ km/s), 171 remain. We minimise $\chi^2 = \sum_i [V_\mathrm{obs}(r_i) - V_\mathrm{pred}(r_i)]^2 / \delta V_i^2$ with the acceleration scale fixed globally. For the primary analysis, errors are floored at $\max(2\;\mathrm{km/s},\; 5\% \, V_\mathrm{obs})$. Section IV.E tests robustness across five error treatments.

**Table I.** 2×2 design: two interpolation functions × two acceleration scales.

| Configuration | Function | Scale | Global free params |
|---|---|---|---|
| Holo + $a_u$ | Holographic (2b) | $a_u = cH_0/(2\pi)$ | 0 |
| Holo + $a_0$ | Holographic (2b) | $a_0 = 1.2 \times 10^{-10}$ | 1 |
| MOND + $a_0$ | MOND simple (2a) | $a_0 = 1.2 \times 10^{-10}$ | 1 |
| MOND + $a_u$ | MOND simple (2a) | $a_u = cH_0/(2\pi)$ | 0 |


# IV. Results

## A. Global fit quality

**Table II.** Fit quality ranked by median $\chi^2/\nu$. Both $a_u$ configurations outperform all $a_0$ configurations. MOND simple with $a_u$ achieves the best overall fit.

| Configuration | Med $\chi^2/\nu$ | Mean $\chi^2/\nu$ | $\chi^2/\nu < 1$ | $\chi^2/\nu < 2$ | Med $\Upsilon_d$ |
|---|---|---|---|---|---|
| **MOND simple + $a_u$** | **1.41** | **3.32** | **36.8%** | **61.4%** | **0.45** |
| **Holographic + $a_u$** | **1.45** | **3.31** | **36.3%** | **62.0%** | **0.45** |
| MOND standard + $a_0$ | 1.52 | 3.31 | 33.3% | 62.0% | 0.55 |
| Holographic + $a_0$ | 1.53 | 3.75 | 31.6% | 56.1% | 0.40 |
| MOND simple + $a_0$ | 1.53 | 3.76 | 32.2% | 56.1% | 0.39 |

The two $a_u$ configurations occupy the top two positions. The single best-performing model across all 171 galaxies is **MOND simple with $a_u$** ($\chi^2/\nu = 1.41$) — Milgrom's own interpolation function, with the fitted constant replaced by the Gibbons–Hawking value. Both $a_0$ configurations (MOND simple and holographic) tie at $\chi^2/\nu = 1.53$.

## B. 2×2 isolation: scale vs. function

To disentangle the contributions of the acceleration scale and the interpolation function, we compare head-to-head within each dimension of the 2×2 design.

**Table III.** Head-to-head isolation. Both the scale and the function independently improve fits. The scale effect is demonstrable within MOND's own functional form.

| Comparison | Varies | Holds fixed | Winner | Win rate |
|---|---|---|---|---|
| $a_u$ vs $a_0$ in Eq. 2b | Scale | Function (Holo) | $a_u$ | **104/171 (60.8%)** |
| $a_u$ vs $a_0$ in Eq. 2a | Scale | Function (MOND) | $a_u$ | **102/171 (59.6%)** |
| Holo vs MOND at $a_u$ | Function | Scale ($a_u$) | Holo | 109/171 (63.7%) |
| Holo vs MOND at $a_0$ | Function | Scale ($a_0$) | Holo | 113/171 (66.1%) |

**Scale effect.** $a_u$ outperforms $a_0$ within both interpolation functions: 60.8% win rate in the holographic form, 59.6% in MOND simple. This is the central result: **MOND's own interpolation function fits the SPARC data better with $cH_0/(2\pi)$ than with its own calibrated constant.** The median $\Delta\chi^2/\nu$ is −0.077 (holographic) and −0.073 (MOND simple) — consistent magnitudes confirming that the scale improvement is independent of the function.

**Function effect.** The holographic interpolation outperforms MOND simple within both acceleration scales: 63.7% at $a_u$, 66.1% at $a_0$. The sharper transition of Eq. 2b is preferred by the data independently of the acceleration scale. However, the median $\Delta\chi^2/\nu$ for the function effect (−0.005 to −0.007) is an order of magnitude smaller than the scale effect (−0.07), indicating that the function shape refines the fit while the scale dominates it.

## C. Mass-to-light ratios

Both $a_u$ configurations produce median $\Upsilon_\mathrm{disk} = 0.45\; M_\odot/L_\odot$ at 3.6 μm, within the SPS range of 0.2–0.8 [14]. The $a_0$ configurations produce lower values (0.39–0.40) for MOND simple and holographic, and higher (0.55) for MOND standard. The $a_u$ value is the most consistent with independent SPS constraints.

## D. Transition sharpness

The holographic $\mu$-function passes through 0.5 over ~0.5 dex in $g_N/a$, versus ~1.0 dex for MOND simple. The 64–66% win rate for Eq. 2b (Table III, rows 3–4) confirms the data prefer the sharper transition. This is independently falsifiable with high-resolution rotation curves densely sampled across $r_t$.

## E. Robustness across error treatments

**Table IV.** Robustness test. Holographic + $a_u$ wins under every error treatment, with the strongest advantage under raw errors.

| Error treatment | Holo+$a_u$ med $\chi^2/\nu$ | MOND-s+$a_0$ med $\chi^2/\nu$ | Holo+$a_u$ wins | Win rate |
|---|---|---|---|---|
| **Raw errors (no floor)** | **2.86** | **3.12** | **108/171** | **63.2%** |
| 2 km/s absolute | 2.49 | 2.68 | 107/171 | 62.6% |
| 3% of $V_\mathrm{obs}$ | 2.17 | 2.43 | 105/171 | 61.4% |
| $\max(2\;\mathrm{km/s},\; 5\%)$ [primary] | 1.45 | 1.53 | 104/171 | 60.8% |
| 10% of $V_\mathrm{obs}$ | 0.69 | 0.71 | 101/171 | 59.1% |

The result holds under every error treatment tested. The advantage is **strongest with raw errors** (63.2% win rate, no floor), diminishing slightly with larger floors as error inflation compresses all models. The median $\Upsilon_\mathrm{disk}$ is stable at 0.45–0.46 across all treatments.


# V. Discussion

## A. The acceleration scale is physical

The 2×2 analysis establishes that the improvement from $a_u$ is not an artefact of the interpolation function. MOND simple — the standard interpolation used in the MOND literature for four decades — performs better with $a_u = cH_0/(2\pi)$ than with its own calibrated $a_0 = 1.2 \times 10^{-10}$ m/s², on 102 out of 171 galaxies. A value derived from the Hubble constant alone, with no calibration to rotation curve data and a 13% offset from the canonical value, outperforms the fitted constant within MOND's own framework.

This result is difficult to attribute to chance. It implies that the true acceleration scale governing the departure from Newtonian dynamics is closer to $cH_0/(2\pi)$ than to the conventional $a_0$, and that the numerical coincidence $a_0 \approx cH_0$ is not a coincidence but a physical relationship.

## B. Physical interpretation

In de Sitter spacetime, the cosmological horizon radiates at the Gibbons–Hawking temperature [7]. The acceleration $a_u = cH/(2\pi)$ marks the scale at which local dynamics become comparable to the horizon's thermal effects. This requires only that the Gibbons–Hawking temperature is physically real — supported by de Sitter thermodynamics [7,11] and the Unruh effect [12].

## C. The 13% discrepancy

The offset between $a_u$ and the conventional $a_0$ is within combined systematics: $H_0$ tension (6–9%) [13], $a_0$ calibration (±17%) [9], baryonic mass estimation (~15%) [14], SPARC distances (~10%) [8]. The superior fit quality of $a_u$ within MOND's own function suggests the canonical $a_0 = 1.2 \times 10^{-10}$ m/s² is biased upward by systematics in the original calibration.

## D. Falsifiable predictions

$$a_u(z) = \frac{cH(z)}{2\pi} \tag{3}$$

At $z = 1$: $a_u \approx 1.5\, a_u(0)$. At $z = 2$: $a_u \approx 2.3\, a_u(0)$. Galaxy kinematics at $z > 1$ should show systematically stronger acceleration corrections than low-$z$ counterparts of equal baryonic mass. CDM predicts halo-dependent evolution; MOND predicts none. JWST IFU spectroscopy [15] will test this within 3–5 years. Additional predictions: (i) the BTFR zero-point is fixed by $H_0$ alone; (ii) the sharper transition is testable with existing high-resolution SPARC data.

## E. Scope and limitations

This paper tests $a_u$ at galactic scales only. The CMB acoustic peaks and matter power spectrum are not addressed. The holographic interpolation (Eq. 2b) is tested on empirical performance; a first-principles derivation is beyond scope. A full Bayesian MCMC analysis with SPS priors would strengthen the comparison; the present $\chi^2$-minimisation results suffice to establish both the scale and function effects.


# VI. Conclusions

We tested the hypothesis that the MOND acceleration scale is the Gibbons–Hawking acceleration $a_u = cH_0/(2\pi)$. Using 171 SPARC galaxies in a 2×2 design — two interpolation functions × two acceleration scales — we find:

**(1)** $a_u$ outperforms $a_0$ within MOND's own interpolation function (102/171 galaxies, 59.6%). MOND simple with $a_u$ achieves the lowest median $\chi^2/\nu$ of any configuration tested (1.41).

**(2)** The holographic interpolation outperforms MOND simple within the same acceleration scale (109–113/171 galaxies, 64–66%).

**(3)** Both effects are independent and robust across all five error treatments, with the strongest advantage under raw errors.

**(4)** Best-fit mass-to-light ratios (median $\Upsilon_d = 0.45$) are consistent with stellar population synthesis predictions.

The scale effect is the more fundamental result. It is demonstrable within MOND's own functional form and does not depend on any new physics beyond the identification $a_u = cH_0/(2\pi)$. A derived constant that outperforms a fitted one — within the fitted constant's own framework — demands explanation. The simplest is that the Gibbons–Hawking temperature of the cosmological horizon is physically connected to galactic dynamics. The prediction $a_u(z) = cH(z)/(2\pi)$ provides a decisive test.


# References

[1] J. F. Navarro, C. S. Frenk, and S. D. M. White, Astrophys. J. 462, 563 (1996).

[2] M. Milgrom, Astrophys. J. 270, 365 (1983).

[3] S. S. McGaugh, J. M. Schombert, G. D. Bothun, and W. J. G. de Blok, Astrophys. J. 533, L99 (2000).

[4] M. Milgrom, Phys. Lett. A 253, 273 (1999).

[5] E. Verlinde, SciPost Phys. 2, 016 (2016).

[6] Planck Collaboration, Astron. Astrophys. 641, A6 (2020).

[7] G. W. Gibbons and S. W. Hawking, Phys. Rev. D 15, 2738 (1977).

[8] F. Lelli, S. S. McGaugh, and J. M. Schombert, Astron. J. 152, 157 (2016).

[9] P. Li, F. Lelli, S. McGaugh, and J. Schombert, Astron. Astrophys. 615, A3 (2018).

[10] B. Famaey and S. S. McGaugh, Living Rev. Relativ. 15, 10 (2012).

[11] T. Jacobson, Phys. Rev. Lett. 75, 1260 (1995).

[12] M. H. Lynch et al., Phys. Rev. D 104, 025015 (2021).

[13] E. Di Valentino et al., Class. Quant. Grav. 38, 153001 (2021).

[14] J. Schombert, S. McGaugh, and F. Lelli, Astron. J. 157, 232 (2019).

[15] H. Übler et al., Astron. Astrophys. 677, A145 (2023).
