# SeMIS: Sequential Multiple Importance Sampling for Reliability Estimation

MATLAB implementation of **Sequential Multiple Importance Sampling (SeMIS)** for estimating small failure probabilities, including in high-dimensional problems. This code accompanies the paper:

> Xia, W., Liao, Z., & Li, B. (2026). Sequential multiple importance sampling for robust and efficient (possibly high-dimensional) reliability estimation. *Mechanical Systems and Signal Processing, 259*, Article 114798. [https://doi.org/10.1016/j.ymssp.2026.114798](https://doi.org/10.1016/j.ymssp.2026.114798)

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

```text
.
├── drive_SeMIS.m            % main driver: selected example, 100 independent runs
├── SeMCS.m                  % sequential sampling loop (Algorithm 1)
├── SeMIS_low_dimension.m    % intermediate distributions, Kriging-based truncation (Eqn. 25)
├── SeMIS_high_dimension.m   % intermediate distributions, PLS-based truncation (Eqn. 26)
├── Support-Fun/MIS.m        % MIS estimator and optional single-run CoV bounds (Eqn. 34–58)
├── SuS.m                    % subset simulation, used as a comparison baseline
├── Limit-state-functions/   % built-in performance functions Ex1 ... Ex12
├── MCMC-samplers/           % aCS, CW_aMH, ES and HMC samplers
├── Support-Fun/             % Nataf transformation, MIS weights and numerical utilities
├── Operations on Logarithmic values/  % log-domain arithmetic
├── LICENSE
└── README.md
```

## Quick start

Edit `example`, `ndim`, `N`, `Nrun` and `SamplerClass` inside `drive_SeMIS.m`. Select `SamplerClass = @aCS` (the default), or choose another supplied sampler: `CW_aMH`, `aCS`, `HMC` or `ES`.

The script defaults to Example 1, 1000 samples per level and 100 runs. `ndim` applies to Examples 10–12; the other examples have fixed dimensions. For a single run, set `Nrun = 1` inside the script.

Run the chosen script from the repository root:

```matlab
drive_SeMIS
```

Use multiple runs to estimate statistics across independent runs.

## Built-in examples

Select an example by setting `example` in `drive_SeMIS.m` to a number from 1 to 12. The examples are numbered as in Tables 1 and 4 of the paper.

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

The script selects the intermediate model from the problem dimension: Kriging below 20 dimensions and PLS otherwise.

## Using a custom performance function

1. Copy one of `Limit-state-functions/Ex*.m` (each is a MATLAB class with the same interface) and rename it, e.g. `MyProblem.m`.
2. Define:
   - `Dist`: probabilistic model of the input vector `X` (see the built-in examples for Nataf-distributed marginals and correlations);
   - `Ndim`, `Nfun`: input dimension and number of performance functions;
   - `g = EvlLSF(obj, u)`: evaluates the performance function on standard-normal samples `u`, where `g < 0` is the failure domain.
3. In `drive_SeMIS.m`, replace the `PF = constructor(...)` assignment inside the run loop with `PF = MyProblem;`. The script selects the intermediate model from `PF.Ndim`.

## Outputs

The raw sampling history `g` contains `Nite` sampled states followed by one final candidate state. `Nsze` records the chain and step counts at each sampled level.

The script leaves `results` in the workspace and prints the mean evaluation count.

| Field in `results` | Meaning |
|---|---|
| `mean`, `standard_deviation` | Per-function mean and sample standard deviation of the failure-probability estimate |
| `cov` | Empirical standard deviation divided by the mean |
| `mean_joint_calls` | Mean number of performance-function evaluations |
| `unfinished` | Number of runs that did not converge |
| `records` | Per-run probabilities, evaluation counts, sample counts and completion flags |


## Optional single-run CoV bounds

The paper (Eqn. 34–58, Appendices A–C) derives analytical bounds of the CoV from a single run, using a three-term decomposition (MCS noise, normalizing-constant noise and cross-term noise). This is disabled by default. To enable it, add the following call after `mis = MIS(out)` inside the run loop:

```matlab
[cv_lw, cv_up] = mis.estimateCov;   % single-run CoV bounds (Eqn. 34–58)
```

## Code-to-paper mapping

| Step | Code | Paper |
|---|---|---|
| Sequence of intermediate distributions | `SeMCS.RunIte` / `UpdObj` | Eqn. (11)–(21), (24)–(27) |
| pCN-MCMC sampling with adaptive step | `MCMC-samplers/aCS.m` | Eqn. (22)–(23), Algorithm 2 |
| Normalizing constants of the ISDs | `MIS.normalizeY` | Eqn. (28)–(32) |
| Balance-heuristic MIS estimator | `MIS` / `WgtMis_Blc` / `MisInt` | Eqn. (4)–(10), (33) |
| CoV bounds (optional) | `MIS.estimateCov` | Eqn. (34)–(58), App. A–C |

## Citation

If you use this toolbox in your research, please cite:

> Xia, W., Liao, Z., & Li, B. (2026). Sequential multiple importance sampling for robust and efficient (possibly high-dimensional) reliability estimation. *Mechanical Systems and Signal Processing, 259*, Article 114798. [https://doi.org/10.1016/j.ymssp.2026.114798](https://doi.org/10.1016/j.ymssp.2026.114798)

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

