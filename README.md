# SeMIS reliability estimation — MATLAB

## Overview
Sequential multiple importance sampling (SeMIS) estimates small failure probabilities using adaptive intermediate distributions and a balance-heuristic MIS estimator. A separate SuS driver uses the classical product estimator. The default sampler is `aCS`; the four available kernels are a unified implementation and do not imply that every kernel/example combination was reported in the cited paper.

## Structure
- `drive_SeMIS.m` / `drive_SuS.m` — independent-run example drivers.
- `./SeMCS.m` — iteration, distribution updates and seed selection.
- `./SeMIS_low_dimension.m` / `SeMIS_high_dimension.m` — GP/PLS-based intermediate distributions and convergence; not iteration drivers.
- `./SuS.m` — SuS intermediate distribution; `Support-Fun/MIS.m` — standalone MIS postprocessor.
- `Limit-state-functions/` — bundled Ex1–Ex12 models.
- `Support-Fun/` — `CW_aMH`, `aCS`, `HMC`, `ES`, shared sampler contract and numerical helpers.

Shared components are vendored inside this repository: no other toolbox checkout is required. Retain the included helper files and license when copying the toolbox.

## Quick start
Run from this repository's root. The explicit small runs below are smoke examples, not the paper's full repetition protocol.

```matlab
results = drive_SeMIS(10, 200, 2, 2026, 20);
baseline = drive_SuS(10, 200, 2, 2026, 20);
```

The driver returns or prints failure-probability estimates and actual model-call counts. Returned `records` include per-run `completed` flags, with an aggregate `unfinished` count. Inspect these flags before treating a run as completed; an estimate alone does not prove convergence.

## Examples and samplers
| Example | Model / dimension | Driver selector |
|---|---|---|
| Ex1 | Piecewise-linear / 2D | `1` |
| Ex2 | Changing topology / 2D | `2` |
| Ex3 | Linear-log / 2D | `3` |
| Ex4 | Polynomial benchmark / 3D | `4` |
| Ex5 | Parallel system / 5D | `5` |
| Ex6 | Brotone Bridge tower / 5D | `6` |
| Ex7 | Noisy performance function / 6D | `7` |
| Ex8 | Reinforced-concrete beam / 7D | `8` |
| Ex9 | Two-DOF damped oscillator / 8D | `9` |
| Ex10 | Linear standard-normal model / configurable n | `10` |
| Ex11 | Nonlinear standard-normal model / configurable n | `11` |
| Ex12 | Linear lognormal model / configurable n | `12` |

- Available samplers: `CW_aMH`, `aCS`, `HMC`, `ES`; default `aCS`.
- Restored defaults: Ex1, 1,000 samples per level, 100 independent repetitions, ambient RNG (no automatic reseeding), conditional probability 0.1, and dimension 2 unless explicitly supplied for Ex10–Ex12. An explicit seed enables deterministic per-run reseeding; otherwise the returned seed field is empty.
- Ex1–Ex9 retain their fixed dimensions. The shared sampler contract uses the model-defined dimension rather than a hard-coded 2D state. Ex10–Ex12 accept the dimension argument; the driver selects GP below 20 dimensions and PLS otherwise. With dimension omitted they therefore use the 2D GP path, not a 50D PLS case. Pass `20`, `50`, `100`, or `1000` explicitly when selecting those high-dimensional paper conditions.
- Select a kernel explicitly:

```matlab
results = drive_SeMIS(10, 200, 2, 2026, 20, 'ES');
```

## Dependencies
MATLAB R2020b or later with Statistics and Machine Learning Toolbox (`fitrgp`, `plsregress`, `makedist`). Supplied drivers are serial and do not require Parallel Computing Toolbox. Ambient serial RNG is not equivalent to the original `parfor` worker streams; native parallel execution has not been validated. Selected high-dimensional and SuS checks used GNU Octave with private compatibility helpers; the public toolbox targets MATLAB. These checks do not establish native MATLAB execution or native `fitrgp` equivalence.

## Method & validation notes
- **Implementation:** SeMIS preserves its existing GP/PLS screening and MIS estimator. SuS retains standard-normal coordinates, Hazen quantiles and thresholds clipped at zero. Failure is strictly `g < 0`. Different toolboxes may select different SuS policies explicitly; common code does not silently unify them.
- **Sampler policy:** aCS starts with unit base sigma and carries its adapted `lambda` in `G` into the next level. Componentwise ES uses an initial angle in `[0, pi]`. HMC uses exact Gaussian-reference Hamiltonian proposals with a likelihood acceptance correction, not generic physical-space gradient HMC. Nataf maps evaluations to physical variables. These are current unified kernels, not a claim of untouched original-paper sampling code.
- **Outputs:** Empirical SD uses independent repeats and sample normalization; CoV is SD divided by the mean. Model-call counts count actual evaluated sample points, not nominal sample totals or screened-out proposals. Optional single-run CoV bounds are a different quantity from empirical repeat-to-repeat CoV.
- **Nataf compatibility:** The author implementation's positive-interval inversion is retained for nonstandard marginal pairs, including its limitations for negative correlations. Normal/lognormal direct transformations and the author's GEV CDF clipping are retained. This release does not silently extend that numerical method.
- **Validation scope:** Focused current-to-cleanup source/trajectory and standalone checks are separate from original-paper statistical comparisons. Native MATLAB low-dimensional `fitrgp` execution was unavailable; no approximate Octave GP was substituted. Interface/trajectory/accounting checks do not establish exact reproduction of every published statistic.
- **Historical scope:** Earlier framework cleanup separated SeMCS, distribution models and MIS. The user's later MATLAB source also changed optional single-run CoV bounds from stochastic normalizer perturbations to analytic bounds. The later author MATLAB implementation is the non-MCMC baseline, not the earliest stochastic-bound version. The recent imposed seed default is reverted; the dimension default is aligned to the author script's value 2. An unqualified claim that all changes since the earliest upload were MCMC-only would be incorrect.
- **Limitations:** Ex8's public reference remains `3.20e-3`, whereas the inspected author-proof table reports `3.01e-3`; neither the model nor its reference was changed to force agreement. Ex10's public `3*sqrt(n)-sum(X)` is a positive scaling of the proof's `3-sum(X)/sqrt(n)`: the zero-threshold failure event agrees, not necessarily adaptive trajectories or costs. The author proof is not a verified final publisher version.

## Reference & license
Weili Xia, Zihan Liao, and Binbin Li. **Sequential multiple importance sampling for robust and efficient (possibly high-dimensional) reliability estimation.** *Mechanical Systems and Signal Processing*, 259 (2026), 114798. [doi:10.1016/j.ymssp.2026.114798](https://doi.org/10.1016/j.ymssp.2026.114798).

Licensed under the [MIT License](LICENSE). The original copyright notice is retained in `LICENSE`.
