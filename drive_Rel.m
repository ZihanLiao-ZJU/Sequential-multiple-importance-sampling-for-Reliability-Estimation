clear;
% -------- Dimension setting (for nD cases only) -------------
Ndim = 20;
% -------- Low dimensional problems -------------------------------
% Tfun = LSF_2D_PiecewiseLinear();                 Pf_ref = 3.20e-5;  % P1
% Tfun = LSF_2D_ChangingToplogical();              Pf_ref = 1.13e-5;  % P2
% Tfun = LSF_2D_LinearLog();                       Pf_ref = 3.23e-5;  % P3
% Tfun = LSF_3D_Benchmarkp6();                     Pf_ref = 3.13e-4;  % P4
% Tfun = LSF_5D_ParallelSystem;                    Pf_ref = 2.13e-4;  % P5
% Tfun = LSF_5D_BrotoneBridge;                     Pf_ref = 2.43e-4;  % P6
% Tfun = LSF_6D_NoisyLSF;                          Pf_ref = 5.29e-4;  % P7
% Tfun = LSF_7D_Reinforced_concrete_beam;          Pf_ref = 3.20e-3;  % P8
% Tfun = LSF_8D_2DOFDampedOcillator;               Pf_ref = 3.78e-7;  % P9

% -------- High dimensional problems ------------------------------
Tfun = LSF_nD_LinearLSF1(Ndim);                  Pf_ref = 1.35e-3;  % P10
% Tfun = LSF_nD_NonlinearLSF1(Ndim);               Pf_ref = 4.75e-6;  % P11
% Tfun = LSF_nD_LinearLSF2(Ndim);                  Pf_ref = 1.91e-3;  % P12

% -------- Engineering problems -------------------------
% Tfun = LSF_10D_Truss;                            Pf_ref = 3.45e-5;  % P13
% Tfun = LSF4_225D_TV_ShearModel;                  Pf_ref = 2.36e-4;  % P14

Ndim = Tfun.Ndim;

% select and initialize intermediate likelihood function and prior distribution
% ----------------------------------------------------------------------------------
% IntBay = SuS(Tfun);
IntBay = SeMIS_low_dimension(Tfun);
% ----------------------------------------------------------------------------------

% put into the constructed sampler
sampler = aCS;

semcs = SeMIS(IntBay,sampler);
semcs.NumSam = 1000;

%% run SuS
Nfun = Tfun.Nfun;
Nrun = 100;

Pf_hat = zeros(Nfun,Nrun);
Ncal   = zeros(1,Nrun);

cv_Pf_lw_all = zeros(Nfun,Nrun);
cv_Pf_up_all = zeros(Nfun,Nrun);

parfor irun = 1:Nrun
    fprintf([num2str(irun),'\n'])
    out_semcs = semcs.RunIte;
    Ncal(irun) = sum(out_semcs.Ncal);
    mis = MIS_Rel(out_semcs);
    % [cv_Pf_lw, cv_Pf_up] = mis.estimateCov;
    % cv_Pf_lw_all(:,irun) = cv_Pf_lw;
    % cv_Pf_up_all(:,irun) = cv_Pf_up;
    Pf_hat(:,irun) = mis.Pf;

    % Pf_hat(:,irun) = IntBay.EvlPf(out_semcs); % SuS
end

Pf_mean = mean(Pf_hat,2);
Pf_std  = std(Pf_hat,0,2);
cv_emp = Pf_std ./ Pf_mean;
cv_lw_mean = mean(cv_Pf_lw_all,2);
cv_up_mean = mean(cv_Pf_up_all,2);

Pf_mu_ref  = mean(Pf_hat./Pf_ref,2);
Pf_sig_ref = std(Pf_hat./Pf_ref,0,2);
Ncal_mu = mean(Ncal);