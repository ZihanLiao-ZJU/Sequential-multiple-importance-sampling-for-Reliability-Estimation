# SeMIS — Sequential Multiple Importance Sampling for Reliability Estimation

MATLAB implementation of the paper

> **Sequential multiple importance sampling for robust and efficient (possibly high-dimensional) reliability estimation**
> (Mechanical Systems and Signal Processing, 2025 — doi: 10.1016/j.ymssp.2024.111950)

The toolbox estimates the failure probability

$$p_f = P\{g(\mathbf X) < 0\} = \int \phi_n(\mathbf u)\, I_{\mathcal F}(\mathbf u)\, \mathrm d\mathbf u$$

by **Sequential Multiple Importance Sampling (SeMIS)**: a sequence of intermediate
distributions is constructed adaptively (Kriging for low-dimensional problems,
PLS for high-dimensional problems), samples are drawn from each with a pMCMC
sampler, and the failure probability is assembled by the balance-heuristic MIS
estimator. A single-run uncertainty quantification (lower/upper bounds of the
coefficient of variation) is available as an option.

---

## Repository layout

```
src v1/
├── drive_SeMIS.m            % main driver: custom example, 100 parallel runs
├── SeMIS.m                  % sequential sampling framework (Algorithm 1)
├── SeMIS_low_dimension.m    % intermediate Bayesian model, Kriging truncation (Eqn. 25)
├── SeMIS_high_dimension.m   % intermediate Bayesian model, PLS truncation (Eqn. 26)
├── MIS_Rel.m                % MIS estimator + optional single-run CoV bounds (Eqn. 34-58)
├── SuS.m                    % subset-simulation interface (comparison baseline)
├── Limit-state-functions/   % built-in performance functions (Ex1 ... Ex12)
└── Support-Fun/             % core support functions (sampler, log-domain ops, ...)
```

## Requirements

- MATLAB R2020b or later
- Statistics and Machine Learning Toolbox (Gaussian-process regression `fitrgp`, distributions)
- Parallel Computing Toolbox (optional; `parfor` degrades to serial without it)

## Quick start

```matlab
>> drive_SeMIS
```

runs Example 1 (`Ex1`) with 100 independent parallel runs and prints the
estimated failure probability, its empirical c.o.v., and the average number of
performance-function evaluations. Edit the `Tfun = Ex1;` line to select any
built-in example (see the list inside `drive_SeMIS.m`).

## Using a custom performance function

1. Copy one of `Limit-state-functions/Ex*.m` (each is a MATLAB class with a
   uniform interface) and rename it, e.g. `MyProblem.m`.
2. Provide:
   - `Dist` — probabilistic model of the input vector `X` (Nataf-distributed marginals/correlations; see the built-in examples),
   - `Ndim`, `Nfun` — input dimension and number of performance functions,
   - `g = EvlLSF(obj, u)` — evaluate the performance function(s) on standard-normal samples `u` (`g < 0` is the failure domain).
3. In `drive_SeMIS.m`, set `Tfun = MyProblem;` (and `Pf_ref` if known), and choose
   `IntBay = SeMIS_low_dimension(Tfun)` (n < 20) or `SeMIS_high_dimension(Tfun)` (n ≥ 20).

## Outputs

With default settings the driver prints, per run and averaged over runs:

- `P_hat_f` — estimated failure probability (Eqn. 9–10, 33)
- `P_hat_f/Pf_ref` — normalized estimate (cf. Tables 2/5 of the paper)
- `CoV (empirical)` — sampled c.o.v. across the independent runs
- `N_cal` — number of performance-function evaluations (tables of the paper)

## Optional single-run CoV bounds

The paper (Eqn. 34–58, Appendices A–C) derives **analytically estimated
bounds of the c.o.v. from a single run** (three-term decomposition: MCS noise +
normalizing-constant noise + cross-term noise). This is disabled by default;
to enable it, uncomment the lines marked `--- CoV estimation (optional) ---`
in `drive_SeMIS.m` (three blocks: allocation, inside the loop, and the print).

```matlab
[cv_lw, cv_up] = mis.estimateCov;   % single-run CoV bounds (Eqn. 34-58)
```

## Method summary (notation of the paper)

| Step | Code | Paper |
|---|---|---|
| Sequence of intermediate distributions | `SeMIS.RunIte` / `UpdObj` | Eqn. (11)–(21), (24)–(27) |
| pMCMC sampling (preconditioned Crank–Nicolson + adaptive step) | `Support-Fun/aCS.m` | Eqn. (22)–(23), Algorithm 2 |
| Normalizing constants of the ISDs | `MIS_Rel.normalizeY` | Eqn. (28)–(32) |
| Balance-heuristic MIS estimator | `MIS_Rel` / `WgtMis_Blc` / `MisInt` | Eqn. (4)–(10), (33) |
| CoV bounds (optional) | `MIS_Rel.estimateCov` | Eqn. (34)–(58), App. A–C |

## Notes

- All probability weights are computed in the log domain
  (`Support-Fun/Operations on Logarithmic values`) for numerical stability.
- Example-specific reference failure probabilities `Pf_ref` (from Tables 1/4 of
  the paper) are listed in the comments of `drive_SeMIS.m`.
- The examples are numbered 1–12 and match Tables 1 and 4 of the paper.

## License

Academic use. Cite the paper above if you use this toolbox in your research.
