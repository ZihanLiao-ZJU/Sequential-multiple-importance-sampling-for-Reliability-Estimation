function results = drive_SuS(example, N, repeats, seed, ndim)
% DRIVE_SUS Run reliability examples through the shared SeMCS engine.
% results = drive_SuS(example, N, repeats, seed, ndim)
% Defaults: Ex1, 1000 samples, 100 runs, seed 2026; ndim=50 for Ex10-12.
% Ex1-9 use their fixed physical dimensions.
if nargin < 1, example = 1; end
if nargin < 2 || isempty(N), N = 1000; end
if nargin < 3, repeats = 100; end
if nargin < 4, seed = 2026; end
if nargin < 5, ndim = 50; end
validateattributes(example, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(repeats, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(seed, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'Limit-state-functions'));
addpath(fullfile(root, 'Support-Fun'));
addpath(fullfile(root, 'Support-Fun', 'Operations on Logarithmic values'));
records = cell(1, repeats);
assert(example <= 12, 'Examples must be in 1:12.');
algorithm = 'SuS';
constructor = str2func(sprintf('Ex%d', example));
for run = 1:repeats
    rng(seed + run - 1, 'twister');
    if example >= 10, PF = constructor(ndim); else, PF = constructor(); end
    model = SuS(PF, 'standard_normal', 'hazen', true);
    sampler = aCS;
    sampler.BayStc = model;
    simulation = SeMCS(sampler);
    simulation.NumSam = N;
    simulation.Verbose = false;
    out = simulation.RunIte;
    probability = 1;
    for level = 1:out.Nite-1
        probability = probability*mean(out.y{level}(4,:) < out.g{level+1}, 'all');
    end
    probability = probability*mean(out.y{end}(4,:) < 0, 'all');
    records{run} = struct('seed', seed + run - 1, ...
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
