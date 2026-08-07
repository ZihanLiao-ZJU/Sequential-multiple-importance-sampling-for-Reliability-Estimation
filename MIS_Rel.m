classdef MIS_Rel
    % MIS_Rel Multiple-importance-sampling reliability estimator (Pf + CoV bounds)
    %
    % Failure-probability estimator (Eqn. 9-10, 33) and single-run uncertainty
    % quantification (Eqn. 34-58 with Appendices A-C):
    %   P_hat_i = (1/N) sum_k W(U_k^(i); Z_hat)                 (Eqn. 33)
    %   Var(P_hat_i) ~= MCS noise + Normalizing constant noise
    %                  + Cross-term noise                         (Eqn. 56)
    %   Var(P_hat_f) in [ sum Var_lower , (sum sqrt(Var_upper))^2 ]  (Eqn. 49/58)

    properties
        TF     % performance function object
        Nite   % number of proposal sets / levels (M+1)
        Pf_hat % failure probability (per PF)  (P_hat_f, Eqn. 9)
        intBay % Bayesian interface providing EvlPDF / EvlLKF
        X
        Y
        G
        log_z  % log normalizing constants: log_z(i) = log z_hat_i (Eqn. 28)
        log_H  % log level-probability estimates: log_H{i} = log H_hat_{i+1} (Eqn. 30/31)

        log_Pi_hat % cell(n_pf,1): store log P_hat_i for each PF (Eqn. 10/33)
    end

    methods
        function obj = MIS_Rel(in)
            obj.Nite   = in.Nite;
            obj.TF     = in.intBay.TF;
            obj.intBay = in.intBay;

            obj.X = in.x;
            obj.Y = in.y;
            obj.G = in.g;

            [obj.log_z, obj.Y, obj.log_H] = obj.normalizeY(obj.X, obj.Y, obj.G);
            log_alpha = WgtMis_Blc(obj.intBay, obj.X, obj.Y, obj.G, obj.log_z);

            n_pf = obj.TF.Nfun;
            Pf_hat = zeros(n_pf, 1);

            % ---- store each log P_hat_i as a property ----
            obj.log_Pi_hat = cell(n_pf, 1);

            for ipf = 1:n_pf
                log_W = cell(obj.Nite, 1);
                for ite = 1:obj.Nite
                    y_i = obj.Y{ite};
                    idx_safe = (y_i(3+ipf, :) >= 0);
                    % log W part before MIS denominator correction (Eqn. 8 numerator)
                    log_W{ite} = y_i(3,:) - y_i(1,:) - y_i(2,:);
                    log_W{ite}(idx_safe) = -inf;      % I_F: g >= 0 -> weight 0
                    log_W{ite} = reshape(log_W{ite}, size(y_i,2:3));
                end

                log_Pi_hat = MisInt(log_alpha, log_W);   % Eqn. 10/33

                % store log(P_hat_i) for this PF
                obj.log_Pi_hat{ipf} = log_Pi_hat;

                Pf_hat(ipf) = exp(logsum(log_Pi_hat, "all"));   % Eqn. 9
            end

            obj.Pf_hat = Pf_hat;
        end

        function [CoV_lw, CoV_up] = estimateCov(obj)
            % Estimate lower and upper bounds of CoV(Pf) from a single run (Eqn. 50-58).
            %
            % The single-run variance is decomposed into the three terms
            % named in Eqn. 55-57:
            %
            %   Var(P_i)_lw = MCS noise + Normalizing constant noise (lower)
            %   Var(P_i)_up = MCS noise + Normalizing constant noise (upper)
            %                 + Cross-term noise
            %
            %   MCS noise                  : Var_U(P_hat_i | Z_hat)          (Eqn. 52/55)
            %   Normalizing constant noise : sum_j sum_l d_ij d_il Cov(Zj,Zl) (Eqn. 51)
            %   Cross-term noise           : sum_j f_ij E[Z_hat_j - z_j]      (Eqn. 55)
            %                                (bias of Z_hat, Eqn. 37 upper)
            %
            % The relative covariance bounds of Z_hat are (Eqn. 38)
            %
            %   Cov_Z_lw(j,l) = sum_{m=1}^{min(j,l)} delta_Hm^2
            %   Cov_Z_up(j,l) = sum_{m=1}^{j} delta_Hm * sum_{t=1}^{l} delta_Ht
            %
            % d_ij and f_ij are estimated by their (C.9) plug-in forms
            %
            %   D_hat_ij = (1/N) sum_k (1/Z_j) W_k alpha_j(U_k),
            %   F_hat_ij = (tau_i/N^2) sum_k 2(W_k - P_i)
            %              * ((1/Z_j) W_k alpha_j,k - D_hat_ij),
            %
            % with alpha_j(U) the balance-heuristic weight (Eqn. 6 / C.3) of
            % the j-th NON-trivial ISD (row j+1 of log_alpha_full; row 1 is
            % the initial level), and Cov(Z_j,Z_l) the absolute covariance
            % zz .* Cov_rel (Eqn. 38).

            n_pf = obj.TF.Nfun;
            M     = obj.Nite - 1;

            % ===============================================================
            % Step 1. Incremental level probabilities H_m and their variances
            % ===============================================================
            H_hat_mean = zeros(M,1);
            var_H  = zeros(M,1);

            for m = 1:M
                H_hat_mean(m) = exp(obj.log_z(m+1) - obj.log_z(m));

                if m == 1
                    h1_hat = H_hat_mean(1);
                    var_H(1) = h1_hat * (1 - h1_hat) / prod(size(obj.Y{1},2:3));
                else
                    [var_H(m), ~] = MoveBlock(obj.log_H{m});
                end
            end

            % delta_Hm^2 = Var(H_hat_m) / H_hat_m^2   (Eqn. 34)
            delta_H_sq = var_H ./ max(H_hat_mean.^2, realmin);
            delta_H_sq = max(delta_H_sq, 0);
            delta_H  = sqrt(delta_H_sq);

            % ===============================================================
            % Step 2. Lower and upper bounds of relative covariance C^Z (Eqn. 38)
            % ===============================================================
            Cov_Z_lw = zeros(M,M);
            Cov_Z_up = zeros(M,M);

            for j = 1:M
                for l = 1:M
                    s_jl = sum(delta_H_sq(1:min(j,l)));

                    Cov_Z_lw(j,l) = s_jl;
                    Cov_Z_up(j,l) = sum(delta_H(1:j)) * sum(delta_H(1:l));
                end
            end

            % ===============================================================
            % Step 3. Normalizing constants
            %
            % z_0 = 1, z_i = exp(log_z(i+1)), i = 1,...,M  (Eqn. 28)
            % ===============================================================
            z_hat   = [1; exp(obj.log_z(2:end))];

            % z_j * z_l factor: absolute covariance Cov(Z_j,Z_l) (Eqn. 38)
            zz = z_hat(2:end) * z_hat(2:end).';

            % ===============================================================
            % Step 4. Bound for |E[Z_hat_j - z_j]|  (Eqn. 37 upper)
            %
            % |E[Z_hat_j - z_j]| <= z_j * sum_{1 <= l < m <= j} delta_Hm delta_Hl
            % ===============================================================
            bias_Z_upper = zeros(M,1);

            for j = 1:M
                delta_j = delta_H(1:j);

                pair_sum = 0.5 * ((sum(delta_j))^2 - sum(delta_j.^2));

                bias_Z_upper(j) = z_hat(j+1) * max(pair_sum, 0);
            end

            % ===============================================================
            % Step 5. MIS denominator part at current run
            % ===============================================================
            [log_alpha, log_alpha_full] = WgtMis_Blc(obj.intBay, obj.X, obj.Y, obj.G, obj.log_z);

            CoV_lw = zeros(n_pf,1);
            CoV_up = zeros(n_pf,1);

            % ===============================================================
            % Step 6. Loop over performance functions
            % ===============================================================
            for ipf = 1:n_pf

                Pf_hat_i = max(obj.Pf_hat(ipf), realmin);

                % (strict (C.9) estimates do not need the component contributions)

                var_mcmc = zeros(obj.Nite,1);
                var_Z_lw = zeros(obj.Nite,1);
                var_Z_up = zeros(obj.Nite,1);
                cross_i  = zeros(obj.Nite,1);   % Cross-term noise (Eqn. 55)

                % -----------------------------------------------------------
                % 6.1 Loop over ISD levels
                % -----------------------------------------------------------
                for iLvl = 1:obj.Nite

                    y_i = obj.Y{iLvl};
                    sz = size(y_i);
                    sz = [sz(2:end), 1];

                    % Failure indicator: I_F = 1 when g(u) < 0
                    idx_safe = (y_i(3 + ipf,:) >= 0);

                    % Basic log weight part before MIS denominator correction (Eqn. 8)
                    log_W_tmp = y_i(3,:) - y_i(1,:) - y_i(2,:);
                    log_W_tmp(idx_safe) = -inf;

                    log_W = zeros(sz);
                    log_W(:) = log_W_tmp(:);

                    % Final empirical MIS log weight:
                    % this is the sequence used in P_hat_i = mean(W_i).
                    log_W_emp = log_alpha{iLvl} + log_W;

                    % -------------------------------------------------------
                    % 6.1.1 MCS noise: Var_U(P_hat_i | Z_hat)  (Eqn. 52/55)
                    %
                    % This follows Eqn. 52 / Appendix B: MoveBlock is applied
                    % to the final empirical MIS weight sequence (the one used
                    % in P_hat_i = mean(W_i)). At the initial level (i = 0)
                    % samples are i.i.d., so Var(W)/N applies (Eqn. B.4).
                    % -------------------------------------------------------
                    if iLvl == 1
                        % i.i.d. samples at the initial level (Eqn. B.4):
                        %   Var_U(P_hat_0 | Z_hat) = Var_U(W) / N
                        W0          = exp(log_W_emp(:));
                        var_mcmc(1) = var(W0) / numel(W0);
                    else
                        [var_mcmc(iLvl), ~] = MoveBlock(log_W_emp);
                    end

                    % Explicit empirical MIS weights for F_hat
                    W     = exp(log_W_emp(:));
                    W_mean = mean(W);

                    % -------------------------------------------------------
                    % 6.1.2 Normalizing constant noise (Eqn. 51)
                    %
                    % Strict (C.9) plug-in estimate of d_ij:
                    %   D_hat_ij = (1/N) sum_k (1/Z_j) W_k alpha_j(U_k),
                    %   alpha_j(U) = balance weight of the j-th non-trivial
                    %   ISD = row j+1 of log_alpha_full. Cov(Z_j,Z_l) is the
                    %   absolute covariance zz .* Cov_rel (Eqn. 38).
                    % -------------------------------------------------------
                    if M > 0
                        D_ij = zeros(M, 1);

                        for j = 1:M
                            alpha_j = exp(log_alpha_full{iLvl}(j+1, :)).';
                            D_ij(j) = mean(W .* alpha_j ./ z_hat(j+1));
                        end

                        var_Z_lw(iLvl) = max(D_ij.' * (zz .* Cov_Z_lw) * D_ij, 0);
                        var_Z_up(iLvl) = max(D_ij.' * (zz .* Cov_Z_up) * D_ij, 0);

                        % ---------------------------------------------------
                        % 6.1.3 Cross-term noise: plug-in F_ij  (Eqn. 54/55, C.9)
                        %
                        % F_ij ~= tau_i/N_i^2 *
                        %        sum 2(W-W_mean)(1/Z_j*W*alpha_j - D_hat)
                        %
                        % MoveBlock estimates Var(P_hat_i | Z_hat), i.e.,
                        %
                        %   Var(P_hat_i | Z_hat)
                        %      ~= (tau_i / N_i^2) * sum_k (W_k - W_mean)^2.
                        %
                        % Hence,
                        %
                        %   tau_i / N_i^2
                        %      ~= Var(P_hat_i | Z_hat) / sum_k (W_k - W_mean)^2.
                        % ---------------------------------------------------
                        ss    = sum((W - W_mean).^2);
                        scale = var_mcmc(iLvl) / max(ss, realmin);

                        F_ij = zeros(M,1);

                        for j = 1:M
                            z_j     = z_hat(j+1);
                            alpha_j = exp(log_alpha_full{iLvl}(j+1, :)).';

                            dW_j = W .* alpha_j ./ z_j;      % (1/Z_j) W alpha_j (C.9)
                            D_hat  = mean(dW_j);             % D_hat_ij (C.9)

                            F_ij(j) = scale * sum( ...
                                2 .* (W - W_mean) .* (dW_j - D_hat) );
                        end

                        cross_i(iLvl) = sum(F_ij .* bias_Z_upper);  % Cross-term noise (Eqn. 55/57, no abs)
                    end
                end

                % -----------------------------------------------------------
                % 6.2 Component variance bounds
                %
                % Var(P_i)_lw = MCS noise + Normalizing constant noise (lower)
                % Var(P_i)_up = MCS noise + Normalizing constant noise (upper)
                %               + Cross-term noise
                % -----------------------------------------------------------
                var_Phat_lw = max(var_mcmc + var_Z_lw, 0);
                var_Phat_up = max(var_mcmc + var_Z_up + cross_i, 0);

                % -----------------------------------------------------------
                % 6.3 Total variance bounds for Pf  (Eqn. 49/58)
                %
                % Lower bound assumes non-negative component cross-covariances.
                % Upper bound follows the Cauchy-Schwarz aggregation.
                % -----------------------------------------------------------
                var_Pf_lw = sum(var_Phat_lw);
                var_Pf_up = (sum(sqrt(var_Phat_up)))^2;

                CoV_lw(ipf) = sqrt(max(var_Pf_lw, 0)) / Pf_hat_i;
                CoV_up(ipf) = sqrt(max(var_Pf_up, 0)) / Pf_hat_i;
            end
        end



        function [log_z, y_norm, log_H] = normalizeY(obj, x, y, g)
            % Sequential estimation of the ISD normalizing constants (Eqn. 28-31):
            %   log_z(i+1) = log_z(i) + logmean( log H_hat_{i+1} )
            % where (nested ISDs) log H_hat_{i+1} = log L_{i+1}(U^(i)) (Eqn. 31),
            % and y_norm stores log L_i shifted by -log z_hat_i (Eqn. 32).
            log_z = zeros(obj.Nite, 1);
            y_norm = y;
            log_H = cell(obj.Nite-1, 1);

            for ite = 1:obj.Nite - 1
                obj.intBay.G = g{ite+1};
                obj.intBay.X = x{ite}(1:end,:);
                obj.intBay.Y = y{ite}(3:end,:);

                log_pi = obj.intBay.EvlPDF;    % log pi_{ite+1}(u)
                log_L  = obj.intBay.EvlLKF;    % log L_{ite+1}(u) (Eqn. 20)

                % log H_hat_{ite+1} = log L_{ite+1} + log pi - log L_ite - log pi_ite
                % (the rPDF pi is unchanged across iterations, so the pi terms cancel)
                log_H{ite} = log_L + log_pi - y{ite}(1, :) - y{ite}(2, :);
                log_H{ite} = log_H{ite} - max(log_H{ite}, [], "all");   % numeric stabilization
                log_z(ite+1) = log_z(ite) + logmean(log_H{ite}, "all");
                y_norm{ite+1}(1,:) = y_norm{ite+1}(1,:) - log_z(ite+1);  % Eqn. 32

                log_H{ite} = reshape(log_H{ite}, size(y{ite},2:3));
                y_norm{ite+1} = reshape(y_norm{ite+1}, size(y{ite+1}));
            end
        end
    end
end
