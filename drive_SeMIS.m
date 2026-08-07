clear;
% =========================================================================
% SeMIS driver
% -------------------------------------------------------------------------
% Runs Sequential Multiple Importance Sampling (SeMIS) on a user-defined
% performance-function example with 100 independent parallel runs and
% estimates the failure probability P_hat_f (Eqn. 9-10, 33).
%
% CoV estimation (single-run uncertainty bounds, Eqn. 34-58) is OPTIONAL
% and disabled by default: uncomment the marked lines below to enable it.
% =========================================================================

% -------- Path setup (run from the src root folder) ---------------------
addpath('Limit-state-functions');
addpath('Support-Fun');
addpath('Support-Fun/Operations on Logarithmic values');

% -------- Custom example ------------------------------------------------
% Low-dimensional examples (paper Table 1, n < 20):
%   Tfun = Ex1;      % piecewise-linear (2D),        Pf_ref = 3.20e-5
%   Tfun = Ex2;      % changing topology (2D),       Pf_ref = 1.13e-5
%   Tfun = Ex3;      % linear-log (2D),              Pf_ref = 3.23e-5
%   Tfun = Ex4;      % polynomial benchmark (3D),    Pf_ref = 3.13e-4
%   Tfun = Ex5;      % parallel system (5D),         Pf_ref = 2.13e-4
%   Tfun = Ex6;      % Brotone Bridge tower (5D),    Pf_ref = 2.370e-4
%   Tfun = Ex7;      % noisy PF (6D),               Pf_ref = 5.29e-4
%   Tfun = Ex8;      % RC beam (7D),                 Pf_ref = 3.20e-3
%   Tfun = Ex9;      % 2-DOF oscillator (8D),        Pf_ref = 3.78e-7
% High-dimensional examples (paper Table 4, n >= 20):
%   Ndim = 50;  Tfun = Ex10(Ndim);   % linear hyperplane,     Pf_ref = 1.35e-3
%   Ndim = 50;  Tfun = Ex11(Ndim);   % nonlinear PF,          Pf_ref = 4.75e-6
%   Ndim = 50;  Tfun = Ex12(Ndim);   % linear (lognormal),     Pf_ref = 1.91e-3
Ndim  = 2;                            % dimension (used by Ex10-12)
Tfun  = Ex1;                          % custom example (default: Example 1)
Pf_ref = 3.20e-5;                     % reference failure probability

% -------- Intermediate Bayesian model (ISD definition) ------------------
% Low-dimensional problems (n < 20):  Kriging-based hybrid truncation (Eqn. 25)
IntBay = SeMIS_low_dimension(Tfun);
% High-dimensional problems (n >= 20): PLS-based hybrid truncation (Eqn. 26)
% IntBay = SeMIS_high_dimension(Tfun);

% -------- Sampler and SeMIS settings ------------------------------------
sampler = aCS;                        % pMCMC sampler (Eqn. 22-23)
semcs = SeMIS(IntBay, sampler);
semcs.NumSam = 1000;                  % sample size per level (N)

% -------- 100 independent parallel runs ---------------------------------
n_pf = Tfun.Nfun;
N_run = 100;

Pf_hat   = zeros(n_pf, N_run);       % P_hat_f per run (Eqn. 9)
N_cal    = zeros(1, N_run);           % PF evaluations per run
% --- CoV estimation (optional; uncomment to enable) ---
% cv_lw_all = zeros(n_pf, N_run);     % estimated CoV lower bounds
% cv_up_all = zeros(n_pf, N_run);     % estimated CoV upper bounds

parfor irun = 1:N_run
    out_semcs = semcs.RunIte;         % sequential sampling (Algorithm 1)
    mis = MIS_Rel(out_semcs);         % MIS estimator (Eqn. 9-10, 33)
    Pf_hat(:, irun)   = mis.Pf_hat;
    % --- CoV estimation (optional; uncomment to enable) ---
    % [cv_lw, cv_up] = mis.estimateCov;  % single-run CoV bounds (Eqn. 34-58)
    % cv_lw_all(:, irun) = cv_lw;
    % cv_up_all(:, irun) = cv_up;
    N_cal(irun) = sum(out_semcs.Ncal);
end

% -------- Output --------------------------------------------------------
Pf_mean  = mean(Pf_hat, 2);           % mean of P_hat_f over runs
Pf_std   = std(Pf_hat, 0, 2);         % std of P_hat_f over runs
cov_emp  = Pf_std ./ Pf_mean;         % empirical (sampled) CoV
N_cal_mean = mean(N_cal);

fprintf('================================================\n');
fprintf('SeMIS results over %d independent runs\n', N_run);
fprintf('  P_hat_f       : %.4e +/- %.4e\n', Pf_mean, Pf_std);
fprintf('  P_hat_f/Pf_ref: %.4f +/- %.4f\n', Pf_mean./Pf_ref, Pf_std./Pf_ref);
fprintf('  CoV (empirical) : %.4f\n', cov_emp);
fprintf('  N_cal (mean)    : %.1f\n', N_cal_mean);
% --- CoV estimation (optional; uncomment to enable) ---
% cv_lw_mean = mean(cv_lw_all, 2);
% cv_up_mean = mean(cv_up_all, 2);
% fprintf('  CoV (estimated) : [%.4f, %.4f]\n', cv_lw_mean, cv_up_mean);
fprintf('================================================\n');
