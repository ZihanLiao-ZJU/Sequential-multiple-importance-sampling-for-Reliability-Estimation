# Shared MATLAB sequential Monte Carlo core

This `semcs-core` branch distributes the shared `+semcs` package under the MIT
license. The twelve MATLAB source files are copied byte for byte from the
validated frozen source-v3 package. Reliability algorithm classes, drivers,
benchmarks and their paper-performance claims are outside this branch.

The package contains one `semcs.Engine`, the model-owned `G` configuration/state
contract, `semcs.MCMCBase`, four Gaussian-reference kernels (`aCS`, `ES`, `HMC`,
`CW_aMH`), the legacy engine adapter, configuration helpers, and the optional
native `semcs.ReSuSProduct` postprocessor including its final true-failure term.
The engine delegates intermediate distributions, thresholds, seed selection and
stopping to the consumer's model. It does not choose a paper's estimator.

Existing consumers resolve an immutable commit through their own `setup_semcs`
and `semcs_dependency.json`. For an explicit local checkout, set `SEMCS_ROOT`
to this directory before calling the consumer's `setup_semcs`. Direct package
users add this directory with `addpath(core_root)` and supply their own model.
No separate reliability-model checkout is required. Algorithm-specific helpers,
targets and native/MIS integrators remain in the consumers.

Models expose `G`, `Ndim`, `Nfun`, `EvlY`, `EvlLKF`, `EvlPDF`, `advance`,
`select_seeds` and `stop`. Optional `SamGen`/`initial_samples` hooks retain native
initial sampling. `G.p`, `G.run`, `G.sampler` and live adapted fields are the
authoritative configuration and state. Engine history preserves packed chain
arrays, actual `Nsam`/`Ncal`, pre-level `g` and `sampled_final_G` snapshots.

MATLAB is required; consumers may additionally require Statistics and Machine
Learning Toolbox or their own model dependencies. HMC here is the exact
Gaussian-reference rotation kernel with likelihood rejection, not a general
gradient-based HMC implementation. Package extraction and same-seed interface
checks do not certify any consumer's agreement with paper statistics.
