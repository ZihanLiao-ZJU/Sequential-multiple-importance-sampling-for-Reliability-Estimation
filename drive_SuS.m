function results = drive_SuS(example, N, repeats, seed, ndim, sampler_name)
% DRIVE_SUS Run reliability examples through the shared SeMCS engine.
% results = drive_SuS(example, N, repeats, seed, ndim, sampler_name)
% sampler_name: CW_aMH, aCS (default), HMC, or ES.
% Defaults: Ex1, 1000 samples, 100 runs, ambient RNG; ndim=2 for Ex10-12.
% Ex1-9 use their fixed physical dimensions.
% Empty/omitted seed preserves ambient RNG; explicit seed resets each run.
% Serial execution does not reproduce original parfor worker random streams.
if nargin < 1, example = 1; end
if nargin < 6 || isempty(sampler_name), sampler_name = 'aCS'; end
assert(any(strcmp(sampler_name,{'CW_aMH','aCS','HMC','ES'})), 'Unknown sampler');
if nargin < 2 || isempty(N), N = 1000; end
if nargin < 3, repeats = 100; end
if nargin < 4, seed = []; end
if nargin < 5, ndim = 2; end
validateattributes(example, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(repeats, {'numeric'}, {'scalar', 'integer', 'positive'});
if ~isempty(seed)
    validateattributes(seed, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
end
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'Limit-state-functions'));
addpath(fullfile(root, 'Support-Fun'));
addpath(fullfile(root, 'Support-Fun', 'Operations on Logarithmic values'));
records = cell(1, repeats);
assert(example <= 12, 'Examples must be in 1:12.');
algorithm = 'SuS';
constructor = str2func(sprintf('Ex%d', example));
for run = 1:repeats
    run_seed = [];
    if ~isempty(seed)
        run_seed = seed + run - 1;
        rng(run_seed, 'twister');
    end
    if example >= 10, PF = constructor(ndim); else, PF = constructor(); end
    model = SuS(PF, 'standard_normal', 'hazen', true);
    sampler = feval(sampler_name);
    sampler.BayStc = model;
    simulation = SeMCS(sampler);
    simulation.NumSam = N;
    simulation.Verbose = false;
    out = simulation.RunIte;
    probability = 1;
    for level = 1:out.Nite-1
        probability = probability*mean(out.y{level}(4,:) < out.g{level+1}.parameters, 'all');
    end
    probability = probability*mean(out.y{end}(4,:) < 0, 'all');
    records{run} = struct('seed', run_seed, ...
        'probability', probability, 'joint_calls', sum(out.Ncal), ...
        'layer_joint_calls', out.Ncal, 'actual_level_samples', out.Nsam, ...
        'completed', out.intBay.FlgCvg);
end
probabilities = cat(2, records{:});
probabilities = cat(2, probabilities.probability);
calls = cellfun(@(record) record.joint_calls, records);
completed = cellfun(@(record) record.completed, records);
sd = std(probabilities, 0, 2);
if repeats < 2, sd(:) = NaN; end
results = struct('example', example, 'algorithm', algorithm, 'N', N, ...
    'repeats', repeats, 'seed', seed, 'mean', mean(probabilities, 2), ...
    'standard_deviation', sd, ...
    'cov', sd ./ mean(probabilities, 2), ...
    'mean_joint_calls', mean(calls), 'unfinished', sum(~completed), ...
    'records', {records});
fprintf('%s example %d: %d runs, mean joint calls %.2f\n', ...
    algorithm, example, repeats, results.mean_joint_calls);
end
