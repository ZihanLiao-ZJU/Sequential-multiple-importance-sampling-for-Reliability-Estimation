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

The common calling sequence is **distribution → sampler → SeMCS → estimator**. `SeMCS` owns the iteration loop; distribution classes define evaluations, intermediate thresholds and convergence; drivers select examples, seed repetitions and summarize results. `MIS` is a separate estimator, not a parent class.

`SuS.m` and `SeMCS.m` are identical in both reliability toolboxes. The SuS drivers choose their numerical policy explicitly: GSS uses physical coordinates, midpoint thresholds and an unclipped terminal threshold; SeMIS uses standard-normal coordinates, Hazen quantiles and clipping at zero. The final failure event is always strictly `g < 0`.

## Requirements

- MATLAB R2020b or later
- Statistics and Machine Learning Toolbox (`fitrgp`, `plsregress`, `makedist`)
- No Parallel Computing Toolbox is required for the supplied drivers.

## Repository structure

```
.
├── drive_SeMIS.m            % main driver: selected example, 100 independent runs
├── SeMCS.m                  % shared sequential sampling framework (Algorithm 1)
├── SeMIS_low_dimension.m    % intermediate distributions, Kriging-based truncation (Eqn. 25)
├── SeMIS_high_dimension.m   % intermediate distributions, PLS-based truncation (Eqn. 26)
├── SuS.m                    % scalar intermediate distribution
├── drive_SuS.m              % independent SuS repetitions
├── Limit-state-functions/   % built-in performance functions Ex1 ... Ex12
├── Support-Fun/             % sampler (aCS), Nataf transformation, MIS estimator and weights, log-domain operations, ...
├── LICENSE
└── README.md
```

## Quick start

Run the driver from the repository root (it adds the subfolders to the MATLAB path):

```matlab
results = drive_SeMIS(1, 1000, 100, 2026);
comparison = drive_SuS(1, 1000, 100, 2026);
```

By default this runs Example 1 (`Ex1`) with 100 independent runs and returns the estimated failure probability, empirical CoV and average model-call count. Use `drive_SeMIS(10, 1000, 100, 2026, 50)` for a 50-dimensional example.

## Built-in examples

Select the example number through the first driver argument. The examples are numbered as in Tables 1 and 4 of the paper.

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

The driver selects `SeMIS_low_dimension` for n < 20 and `SeMIS_high_dimension` otherwise. The fifth argument sets the dimension of `Ex10`–`Ex12`.

## Using a custom performance function

1. Copy one of `Limit-state-functions/Ex*.m` (each is a MATLAB class with the same interface) and rename it, e.g. `MyProblem.m`.
2. Define:
   - `Dist`: probabilistic model of the input vector `X` (see the built-in examples for Nataf-distributed marginals and correlations);
   - `Ndim`, `Nfun`: input dimension and number of performance functions;
   - `g = EvlLSF(obj, x)`: evaluates the performance function on physical samples `x` (the model performs the Nataf transformation), where `g < 0` is the failure domain.
3. Construct `model = SeMIS_low_dimension(MyProblem)` for n < 20 or `model = SeMIS_high_dimension(MyProblem)` otherwise. Set `sampler = aCS; sampler.BayStc = model; simulation = SeMCS(sampler);`, set `simulation.NumSam`, then run `out = simulation.RunIte; estimate = MIS(out);`.

## Outputs

All drivers return the same summary fields: `example`, `algorithm`, `N`, `repeats`, `seed`, `mean`, `standard_deviation`, `cov`, `mean_joint_calls`, `unfinished`, and `records`.

Each record contains `seed`, `probability`, `joint_calls`, `layer_joint_calls`, `actual_level_samples`, and `completed`. `unfinished` is the number of repetitions that did not reach the failure domain. `cov` is the empirical sample standard deviation divided by the mean; it is distinct from optional single-run uncertainty bounds. With one repetition, sample standard deviation and empirical CoV are returned as `NaN`.

The intermediate sampling output from `SeMCS.RunIte` retains `x`, `y`, `g`, `Nsam`, `Nsze`, `Ncal`, `Nite`, `flg_cvg`, and the final distribution object `intBay`. GSS also retains the active response history `h`.

## Optional single-run CoV bounds

The paper derives optional single-run CoV bounds. After obtaining a sampling output, use the separate estimator directly:

```matlab
mis = MIS(out);
[cv_lw, cv_up] = mis.estimateCov;   % single-run CoV bounds (Eqn. 34–58)
```

## Code-to-paper mapping

| Step | Code | Paper |
|---|---|---|
| Sequence of intermediate distributions | `SeMCS.RunIte` / `UpdObj` | Eqn. (11)–(21), (24)–(27) |
| pCN-MCMC sampling with adaptive step | `Support-Fun/aCS.m` | Eqn. (22)–(23), Algorithm 2 |
| Normalizing constants of the ISDs | `MIS.normalizeY` | Eqn. (28)–(32) |
| Balance-heuristic MIS estimator | `MIS` / `WgtMis_Blc` / `MisInt` | Eqn. (4)–(10), (33) |
| CoV bounds (optional) | `MIS.estimateCov` | Eqn. (34)–(58), App. A–C |

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
