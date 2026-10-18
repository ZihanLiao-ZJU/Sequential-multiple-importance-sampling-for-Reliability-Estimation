# SeMIS: Sequential Multiple Importance Sampling for Reliability Estimation

MATLAB implementation of **Sequential Multiple Importance Sampling (SeMIS)** for estimating small failure probabilities, including in high-dimensional problems. This code accompanies the paper:

> **Sequential multiple importance sampling for robust and efficient (possibly high-dimensional) reliability estimation.**
> *Mechanical Systems and Signal Processing*, 259 (2026) 114798. doi: [10.1016/j.ymssp.2026.114798](https://doi.org/10.1016/j.ymssp.2026.114798)

## Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Repository structure](#repository-structure)
- [Quick start](#quick-start)
- [Built-in examples](#built-in-examples)
- [Using a custom performance function](#using-a-custom-performance-function)
- [Outputs](#outputs)
- [Optional single-run CoV bounds](#optional-single-run-cov-bounds)
- [Code-to-paper mapping](#code-to-paper-mapping)
- [Citation](#citation)
- [License](#license)

## Overview

SeMIS estimates the failure probability

$$p_f = P\{g(\mathbf X) < 0\} = \int \phi_n(\mathbf u)\, I_{\mathcal F}(\mathbf u)\, \mathrm d\mathbf u$$

by constructing a sequence of intermediate importance sampling distributions (ISDs) that move progressively toward the failure domain:

1. The intermediate distributions are defined adaptively, using a Kriging-based truncation for low-dimensional problems or a PLS-based truncation for high-dimensional problems.
2. Samples are drawn from each intermediate distribution with a preconditioned Crank–Nicolson (pCN) MCMC sampler with adaptive step size.
3. All samples are combined with the balance-heuristic multiple importance sampling (MIS) estimator to obtain the failure probability.

As an option, lower and upper bounds of the coefficient of variation (CoV) can be estimated from a single run. All probability weights are computed in the log domain for numerical stability.

## Requirements

- MATLAB R2020b or later
- Statistics and Machine Learning Toolbox (`fitrgp`, `plsregress`, `makedist`)
- Parallel Computing Toolbox (optional; without it, `parfor` runs serially)

## Repository structure

```
.
├── drive_SeMIS.m            % main driver: selected example, 100 independent runs
├── SeMIS.m                  % sequential sampling framework (Algorithm 1)
├── SeMIS_low_dimension.m    % intermediate distributions, Kriging-based truncation (Eqn. 25)
├── SeMIS_high_dimension.m   % intermediate distributions, PLS-based truncation (Eqn. 26)
├── MIS_Rel.m                % MIS estimator and optional single-run CoV bounds (Eqn. 34–58)
├── SuS.m                    % subset simulation, used as a comparison baseline
├── Limit-state-functions/   % built-in performance functions Ex1 ... Ex12
├── Support-Fun/             % sampler (aCS), Nataf transformation, MIS weights, log-domain operations, ...
├── LICENSE
└── README.md
```

## Quick start

Run the driver from the repository root (it adds the subfolders to the MATLAB path):

```matlab
>> drive_SeMIS
```

By default this runs Example 1 (`Ex1`) with 100 independent runs and prints the estimated failure probability, its empirical CoV, and the average number of performance-function evaluations.

## Built-in examples

Select an example by editing the `Tfun = ...` line in `drive_SeMIS.m` and set `Pf_ref` accordingly. The examples are numbered as in Tables 1 and 4 of the paper.

| Example | Description | Dimension | `Pf_ref` |
|---|---|---|---|
| `Ex1` | piecewise-linear | 2 | 3.20e-5 |
| `Ex2` | changing topology | 2 | 1.13e-5 |
| `Ex3` | linear-log | 2 | 3.23e-5 |
| `Ex4` | polynomial benchmark | 3 | 3.13e-4 |
| `Ex5` | parallel system | 5 | 2.13e-4 |
| `Ex6` | Brotone Bridge tower | 5 | 2.370e-4 |
| `Ex7` | noisy performance function | 6 | 5.29e-4 |
| `Ex8` | RC beam | 7 | 3.20e-3 |
| `Ex9` | 2-DOF oscillator | 8 | 3.78e-7 |
| `Ex10(50)` | linear hyperplane | 50 | 1.35e-3 |
| `Ex11(50)` | nonlinear performance function | 50 | 4.75e-6 |
| `Ex12(50)` | linear, lognormal inputs | 50 | 1.91e-3 |

For the high-dimensional examples (`Ex10`–`Ex12`), also switch the intermediate model in `drive_SeMIS.m` to `IntBay = SeMIS_high_dimension(Tfun);`.

## Using a custom performance function

1. Copy one of `Limit-state-functions/Ex*.m` (each is a MATLAB class with the same interface) and rename it, e.g. `MyProblem.m`.
2. Define:
   - `Dist`: probabilistic model of the input vector `X` (see the built-in examples for Nataf-distributed marginals and correlations);
   - `Ndim`, `Nfun`: input dimension and number of performance functions;
   - `g = EvlLSF(obj, u)`: evaluates the performance function on standard-normal samples `u`, where `g < 0` is the failure domain.
3. In `drive_SeMIS.m`, set `Tfun = MyProblem;` (and `Pf_ref`, if known), then choose `IntBay = SeMIS_low_dimension(Tfun)` for n < 20 or `IntBay = SeMIS_high_dimension(Tfun)` for n ≥ 20.

## Outputs

The driver prints the following quantities, averaged over the independent runs:

| Output | Meaning |
|---|---|
| `P_hat_f` | estimated failure probability (Eqn. 9–10, 33) |
| `P_hat_f/Pf_ref` | estimate normalized by the reference value (cf. Tables 2 and 5 of the paper) |
| `CoV (empirical)` | coefficient of variation across the independent runs |
| `N_cal (mean)` | average number of performance-function evaluations per run |

## Optional single-run CoV bounds

The paper (Eqn. 34–58, Appendices A–C) derives analytical bounds of the CoV from a single run, using a three-term decomposition (MCS noise, normalizing-constant noise and cross-term noise). This is disabled by default. To enable it, uncomment the three blocks marked `--- CoV estimation (optional) ---` in `drive_SeMIS.m` (allocation, inside the loop, and the print statement):

```matlab
[cv_lw, cv_up] = mis.estimateCov;   % single-run CoV bounds (Eqn. 34–58)
```

## Code-to-paper mapping

| Step | Code | Paper |
|---|---|---|
| Sequence of intermediate distributions | `SeMIS.RunIte` / `UpdObj` | Eqn. (11)–(21), (24)–(27) |
| pCN-MCMC sampling with adaptive step | `Support-Fun/aCS.m` | Eqn. (22)–(23), Algorithm 2 |
| Normalizing constants of the ISDs | `MIS_Rel.normalizeY` | Eqn. (28)–(32) |
| Balance-heuristic MIS estimator | `MIS_Rel` / `WgtMis_Blc` / `MisInt` | Eqn. (4)–(10), (33) |
| CoV bounds (optional) | `MIS_Rel.estimateCov` | Eqn. (34)–(58), App. A–C |

## Citation

If you use this code in your research, please cite:

> Xia, W., Liao, Z., & Li, B. (2026). Sequential multiple importance sampling for robust and efficient (possibly high-dimensional) reliability estimation. *Mechanical Systems and Signal Processing*, 259, 114798. https://doi.org/10.1016/j.ymssp.2026.114798

```bibtex
@article{xia2026sequential,
  title   = {Sequential multiple importance sampling for robust and efficient (possibly high-dimensional) reliability estimation},
  author  = {Xia, Weili and Liao, Zihan and Li, Binbin},
  journal = {Mechanical Systems and Signal Processing},
  volume  = {259},
  pages   = {114798},
  year    = {2026},
  doi     = {10.1016/j.ymssp.2026.114798}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](./LICENSE) for details.
