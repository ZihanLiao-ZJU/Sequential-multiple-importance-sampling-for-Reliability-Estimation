# Sequential Multiple Importance Sampling (SeMIS)

MATLAB research code for:

> W. Xia, Z. Liao, and B. Li, "Sequential multiple importance sampling for robust and efficient (possibly high-dimensional) reliability estimation," *Mechanical Systems and Signal Processing*, Article 114798, 2026. [https://doi.org/10.1016/j.ymssp.2026.114798](https://doi.org/10.1016/j.ymssp.2026.114798)

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
