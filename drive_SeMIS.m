%% Problem and sampling settings
PF = Ex1();
% PF = Ex10(50);             % Examples 10–12 accept a dimension
N = 1000;                       % Constant samples per level
Nrun = 100;                   % Independent runs; set to 1 for a single run
SamplerClass = @aCS;      % CW_aMH, aCS, HMC or ES

%% Model and sampler
if PF.Ndim < 20
    model = SeMIS_low_dimension(PF);
else
    model = SeMIS_high_dimension(PF);
end
simulation = SeMCS(SamplerClass(model));
simulation.NumSam = N;

%% Run simulations
records = cell(1, Nrun);
algorithm = 'SeMIS';
problem_name = class(PF);
example = sscanf(problem_name, 'Ex%d');   % Built-in benchmark number in results
for run = 1:Nrun
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
fprintf('%s %s: %d runs, mean joint calls %.2f\n', ...
    algorithm, problem_name, Nrun, results.mean_joint_calls);
