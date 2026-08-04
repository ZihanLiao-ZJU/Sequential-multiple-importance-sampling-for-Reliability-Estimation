classdef MIS_Rel
    % MIS_Rel Multiple-importance-sampling reliability estimator (Pf + CoV bounds)

    properties
        TF     % limit-state function object
        Nite   % number of proposal sets / levels
        Pf     % failure probability (per LSF)
        intBay % Bayesian interface providing EvlPDF / EvlLKF
        X
        Y
        G
        logC
        W_res

        Pi_log % cell(Nlsf,1): store pi_log for each LSF
    end

    methods
        function obj = MIS_Rel(in)
            obj.Nite   = in.Nite;
            obj.TF     = in.intBay.TF;
            obj.intBay = in.intBay;

            obj.X = in.x;
            obj.Y = in.y;
            obj.G = in.g;

            [obj.logC, obj.Y, obj.W_res] = obj.normalizeY(obj.X, obj.Y, obj.G);
            w_mis = WgtMis_Blc(obj.intBay, obj.X, obj.Y, obj.G, obj.logC);

            Nlsf = obj.TF.Nfun;
            Pf_loc = zeros(Nlsf, 1);

            % ---- store each pi_log as a property ----
            obj.Pi_log = cell(Nlsf, 1);

            for ilsf = 1:Nlsf
                w_is = cell(obj.Nite, 1);
                for ite = 1:obj.Nite
                    yi = obj.Y{ite};
                    idx_safe = (yi(3+ilsf, :) >= 0);
                    w_is{ite} = yi(3,:) - yi(1,:) - yi(2,:);
                    w_is{ite}(idx_safe) = -inf;
                    w_is{ite} = reshape(w_is{ite}, size(yi,2:3));
                end

                pi_log = MisInt(w_mis, w_is);

                % store log(pi) for this LSF
                obj.Pi_log{ilsf} = pi_log;

                Pf_loc(ilsf) = exp(logsum(pi_log, "all"));
            end

            obj.Pf = Pf_loc;
        end

        function [CoV_Pf_lw, CoV_Pf_up] = estimateCov(obj)
            N_lsf = obj.TF.Nfun;
            M     = obj.Nite - 1;
            N_ol  = 1e3;
            % ---------- Var(Hhat_m) ----------
            Var_Hhat = zeros(M,1);
            for m = 1:M
                if m == 1
                    p = exp(obj.logC(2));
                    Var_Hhat(1) = p*(1-p) / prod(size(obj.Y{1},2:3));
                else
                    [Var_Hhat(m),~] = MoveBlock(obj.W_res{m});
                end
            end
            % ---------- mean(Hhat_m) ----------
            Hhat_mean = zeros(M,1);
            for m = 1:M
                Hhat_mean(m) = exp(obj.logC(m+1) - obj.logC(m));
            end
            % ==============================================================
            % LOWER: Corr(Hhat)=0
            % ==============================================================
            Hhat_smp = zeros(M,N_ol);
            for m = 1:M
                Hhat_smp(m,:) = sqrt(Var_Hhat(m)).*randn(1,N_ol) + Hhat_mean(m);
            end
            logZhat_smp = log(cumprod(Hhat_smp,1));
            % ==============================================================
            % UPPER: Corr(Hhat)=1
            % ==============================================================
            z0 = randn(1,N_ol);
            Hhat_smp_up = zeros(M,N_ol);
            for m = 1:M
                Hhat_smp_up(m,:) = sqrt(Var_Hhat(m)).*z0 + Hhat_mean(m);
            end
            Hhat_smp_up = max(Hhat_smp_up, realmin);
            logZhat_smp_up = log(cumprod(Hhat_smp_up,1));
            CoV_Pf_lw = zeros(N_lsf,1);
            CoV_Pf_up = zeros(N_lsf,1);

            for ilsf = 1:N_lsf
                % ================= LOWER =================
                logPhat_i_smp = zeros(obj.Nite, N_ol);
                Var_inner_i_smp = zeros(obj.Nite, N_ol);

                for r = 1:N_ol
                    y_rnd = obj.Y;
                    for m = 1:M
                        y_rnd{m+1}(1,:) = obj.Y{m+1}(1,:) + obj.logC(m+1) - logZhat_smp(m,r);
                    end
                    w_mis = WgtMis_Blc(obj.intBay, obj.X, obj.Y, obj.G, [0; logZhat_smp(:,r)]);
                    w_is = cell(obj.Nite,1);
                    for iLvl = 1:obj.Nite
                        yi = y_rnd{iLvl};
                        sz = size(yi); sz = [sz(2:end), 1];

                        idx_safe = (yi(3+ilsf,:) >= 0);
                        w_tmp = yi(3,:) - yi(1,:) - yi(2,:);
                        w_tmp(idx_safe) = -inf;

                        w_is{iLvl} = zeros(sz);
                        w_is{iLvl}(:) = w_tmp(:);
                    end
                    logPhat_i_smp(:,r) = MisInt(w_mis, w_is);

                    for iLvl = 1:obj.Nite
                        pi_est = w_mis{iLvl} + w_is{iLvl};
                        if iLvl == 1
                            pi_est = reshape(pi_est, [], 10);
                        end
                        [Var_inner_i_smp(iLvl,r),~] = MoveBlock(pi_est);
                    end
                end
                Var_outer_i = var(exp(logPhat_i_smp), 0, 2);
                Var_Phat_i  = Var_outer_i + mean(Var_inner_i_smp,2);
                CoV_Pf_lw(ilsf) = sqrt(sum(Var_Phat_i)) / obj.Pf(ilsf);

                % ================= UPPER =================
                logPhat_i_smp(:) = 0;
                Var_inner_i_smp(:) = 0;

                for r = 1:N_ol
                    y_rnd = obj.Y;
                    for m = 1:M
                        y_rnd{m+1}(1,:) = obj.Y{m+1}(1,:) + obj.logC(m+1) - logZhat_smp_up(m,r);
                    end
                    w_mis = WgtMis_Blc(obj.intBay, obj.X, obj.Y, obj.G, [0; logZhat_smp_up(:,r)]);
                    w_is = cell(obj.Nite,1);
                    for iLvl = 1:obj.Nite
                        yi = y_rnd{iLvl};
                        sz = size(yi); sz = [sz(2:end), 1];
                        idx_safe = (yi(3+ilsf,:) >= 0);
                        w_tmp = yi(3,:) - yi(1,:) - yi(2,:);
                        w_tmp(idx_safe) = -inf;
                        w_is{iLvl} = zeros(sz);
                        w_is{iLvl}(:) = w_tmp(:);
                    end
                    logPhat_i_smp(:,r) = MisInt(w_mis, w_is);
                    for iLvl = 1:obj.Nite
                        pi_est = w_mis{iLvl} + w_is{iLvl};
                        if iLvl == 1
                            pi_est = reshape(pi_est, [], 10);
                        end
                        [Var_inner_i_smp(iLvl,r),~] = MoveBlock(pi_est);
                    end
                end
                Var_outer_i = var(exp(logPhat_i_smp), 0, 2);
                Var_Phat_i  = Var_outer_i + mean(Var_inner_i_smp,2);
                CoV_Pf_up(ilsf) = sqrt((sum(sqrt(Var_Phat_i)))^2) / obj.Pf(ilsf);
            end
        end

        function [c_log, y_nor, w_res] = normalizeY(obj, x, y, g)
            c_log = zeros(obj.Nite, 1);
            y_nor = y;
            w_res = cell(obj.Nite-1, 1);

            for ite = 1:obj.Nite - 1
                obj.intBay.G = g{ite+1};
                obj.intBay.X = x{ite}(1:end,:);
                obj.intBay.Y = y{ite}(3:end,:);

                pi_pdf = obj.intBay.EvlPDF;
                Li_lkf = obj.intBay.EvlLKF;

                w_res{ite} = Li_lkf + pi_pdf - y{ite}(1, :) - y{ite}(2, :);
                w_res{ite} = w_res{ite} - max(w_res{ite}, [], "all");
                c_log(ite+1) = c_log(ite) + logmean(w_res{ite}, "all");
                y_nor{ite+1}(1,:) = y_nor{ite+1}(1,:) - c_log(ite+1);

                w_res{ite} = reshape(w_res{ite}, size(y{ite},2:3));
                y_nor{ite+1} = reshape(y_nor{ite+1}, size(y{ite+1}));
            end
        end
    end
end