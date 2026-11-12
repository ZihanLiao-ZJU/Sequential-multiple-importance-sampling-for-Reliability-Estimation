%% Problem and sampling settings
example = 1;
N = 1000;                       % Constant samples per level
Nrun = 100;                   % Independent runs; set to 1 for a single run
SamplerClass = @aCS;      % CW_aMH, aCS, HMC or ES
ndim = 2;                    % Dimension for Examples 10-12

%% Run simulations
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'Limit-state-functions'));
addpath(fullfile(root, 'Support-Fun'));
addpath(fullfile(root, 'MCMC-samplers'));
addpath(fullfile(root, 'Operations on Logarithmic values'));
records = cell(1, Nrun);
assert(example <= 12, 'Examples must be in 1:12.');
algorithm = 'SeMIS';
constructor = str2func(sprintf('Ex%d', example));
for run = 1:Nrun
    if example >= 10, PF = constructor(ndim); else, PF = constructor(); end
    if PF.Ndim < 20
        model = SeMIS_low_dimension(PF);
    else
        model = SeMIS_high_dimension(PF);
    end
    sampler = SamplerClass();
    sampler.BayStc = model;
    simulation = SeMCS(sampler);
    simulation.NumSam = N;
    simulation.Verbose = false;
    out = simulation.RunIte;
    mis = MIS(out);
    probability = mis.Pf_hat;
    records{run} = struct('probability', probability, 'joint_calls', sum(out.Ncal), ...
        'layer_joint_calls', out.Ncal, 'actual_level_samples', out.Nsam, ...
        'completed', out.intBay.FlgCvg);
end
probabilities = cat(2, records{:});
probabilities = cat(2, probabilities.probability);
calls = cellfun(@(record) record.joint_calls, records);
completed = cellfun(@(record) record.completed, records);
sd = std(probabilities, 0, 2);
if Nrun < 2, sd(:) = NaN; end
results = struct('example', example, 'algorithm', algorithm, 'N', N, ...
    'repeats', Nrun, 'mean', mean(probabilities, 2), ...
    'standard_deviation', sd, ...
    'cov', sd ./ mean(probabilities, 2), ...
    'mean_joint_calls', mean(calls), 'unfinished', sum(~completed), ...
    'records', {records});
fprintf('%s example %d: %d runs, mean joint calls %.2f\n', ...
    algorithm, example, Nrun, results.mean_joint_calls);
