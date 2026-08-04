# Sequential Multiple Importance Sampling (SeMIS)

MATLAB research code for:

> W. Xia, Z. Liao, and B. Li, "Sequential multiple importance sampling for robust and efficient (possibly high-dimensional) reliability estimation," *Mechanical Systems and Signal Processing*, Article 114798, 2026. [https://doi.org/10.1016/j.ymssp.2026.114798](https://doi.org/10.1016/j.ymssp.2026.114798)

> **Status:** The core SeMIS algorithm is implemented, but this repository is not yet a complete one-to-one reproduction package for all results in the paper. See [Known differences](#known-differences).

## Overview

SeMIS estimates small failure probabilities by combining samples from a sequence of intermediate importance sampling densities. The implementation includes:

- isotropic variance inflation for global exploration;
- hybrid truncation to reduce expensive performance-function evaluations;
- Kriging screening for problems with dimension `n < 20`;
- one-component PLS screening for problems with dimension `n >= 20`;
- adaptive conditional sampling using a pCN-type proposal;
- balance-heuristic multiple importance sampling using samples from all iterations;
- Nataf transformation between the physical and standard normal spaces.

## Main files

| File | Purpose |
|---|---|
| `SeMIS.m` | Main sequential workflow |
| `SeMIS_low_dimension.m` | Variance inflation and Kriging-based truncation |
| `SeMIS_high_dimension.m` | Variance inflation and PLS-based truncation |
| `MIS_Rel.m` | MIS failure-probability estimator |
| `SuS.m` | Subset simulation baseline |
| `Support-Fun/aCS.m` | Adaptive pCN conditional sampler |
| `Support-Fun/Nataf.m` | Nataf transformation |
| `Limit-state-functions/` | Benchmark and engineering examples |
| `drive_Rel.m` | Batch experiment script |
| `plot/` | Plotting scripts |

## Requirements

- MATLAB with `pagemtimes` support;
- Statistics and Machine Learning Toolbox (`fitrgp`, `plsregress`, `makedist`, and `randsample`);
- Parallel Computing Toolbox for scripts using `parfor`.

The following dependencies are not included:

- `FEM`, required by Example 13;
- `GauIntPot`, required by the general correlated non-normal branch of `Nataf`;
- the 225-dimensional structural model used in Example 14.

## Quick start

```matlab
clear; clc;
addpath(genpath(pwd));
rng(1);

Tfun = LSF_2D_PiecewiseLinear();

if Tfun.Ndim < 20
    intBay = SeMIS_low_dimension(Tfun);
else
    intBay = SeMIS_high_dimension(Tfun);
end

semis = SeMIS(intBay, aCS);
semis.NumSam = 1000;

out = semis.RunIte;
mis = MIS_Rel(out);

Pf_hat = mis.Pf
Ncal   = sum(out.Ncal)
```

Important output fields are:

- `mis.Pf`: estimated failure probability;
- `out.x`, `out.y`, and `out.g`: samples, evaluations, and intermediate-density parameters;
- `out.Ncal`: performance-function evaluations per iteration;
- `out.Nite`: number of iterations;
- `mis.Pi_log`: contribution of each intermediate density to the MIS estimate.

## Paper settings

The main settings used in the paper are:

- `NumSam = 1000` per iteration;
- intermediate probability `p0 = 0.1`;
- Kriging for `n < 20` and PLS for `n >= 20`;
- 1000 independent runs for Examples 1-12;
- 100 independent runs for Examples 13-14.

## Known differences

- Example 10 uses `beta = 6` in the current code, while the paper uses `beta = 3`.
- `drive_Rel.m` currently selects the low-dimensional implementation at `n = 20`; the paper uses PLS for `n >= 20`.
- `drive_Rel.m` uses 100 runs instead of 1000 for Examples 1-12.
- The current `MIS_Rel.estimateCov` is not a direct implementation of the final Taylor-expansion formulas in Section 4.2 and Appendices A-C.
- Example 8 currently corresponds to a failure probability near `3.20e-3`, while the paper reports `3.01e-3`.
- Some reference probabilities in `drive_Rel.m`, including Examples 6 and 12, require synchronization with the paper.
- Example 13 cannot run without the missing `FEM` class; Example 14 is not included.
- Several plotting scripts contain hard-coded results or data from earlier paper versions.
- The minus sign shown for Example 4 in the uncorrected proof is likely a typesetting error. The plus sign used in the code is consistent with the reported reference probability.

These differences are expected to be updated in later revisions.

## Citation

```bibtex
@article{Xia2026SeMIS,
  title   = {Sequential multiple importance sampling for robust and efficient (possibly high-dimensional) reliability estimation},
  author  = {Xia, Weili and Liao, Zihan and Li, Binbin},
  journal = {Mechanical Systems and Signal Processing},
  year    = {2026},
  pages   = {114798},
  doi     = {10.1016/j.ymssp.2026.114798}
}
```
