classdef aCS
    % Adaptive conditional sampling (pCN-type MCMC) with Nataf transformation
    properties(Constant)
        % SedSav  : use provided seeds at k=1
        % AcpOpt  : target acceptance rate
        % NumIte  : number of adaptive rounds
        % Lambda  : initial proposal scaling
        SedSav = false;
        AcpOpt = 0.44;
        NumIte = 5;
        Lambda = 0.6;
    end

    properties
        % BayStc : Bayesian structure (likelihood, prior, counters, etc.)
        BayStc
    end

    methods
        function [u,y,N_cal,baystc] = SamGen(obj,u_seed,y_seed,Ns)
            % Generate MCMC samples targeting q_i = pi_i * L_i / z_i  (Eqn. 21)
            % pCN proposal (Eqn. 22) + adaptive scaling on acceptance rate.

            % initialization
            baystc = obj.BayStc;
            sedsav = obj.SedSav;
            Ndim = obj.BayStc.Ndim;
            Nfun = obj.BayStc.Nfun;
            n_round = obj.NumIte;
            lambda_pcn = obj.Lambda;
            acp_opt = obj.AcpOpt;

            % Cholesky factor of the rPDF covariance (Eqn. 12: sigma^2 I)
            if isfield(baystc.G,'L_pi')
                L_pi = baystc.G.L_pi;
            else
                L_pi = eye(Ndim);
            end

            Nc = size(u_seed,2);
            N_cal = 0;

            % allocate memory
            u = zeros(Ndim,Nc,Ns);
            y = zeros(Nfun,Nc,Ns);

            % Gaussian innovations (Eqn. 22: v ~ N(0,I)) and MH uniforms
            v = pagemtimes(L_pi,randn(Ndim,Nc,Ns));
            u_mh = rand(1,Nc,Ns);

            % group chains for adaptive updates
            n_grp = max(ceil(Nc/n_round),1);
            n_round = ceil(Nc/n_grp);
            idx_grp = min(repmat([1,n_grp],n_round,1)+(0:n_grp:(n_round-1)*n_grp)',Nc);

            % base proposal scale
            % sigma0 = max(std(u_seed,0,2),eps);
            sigma0 = 1;

            for iter = 1:n_round
                % pCN parameters (Eqn. 22): xi = rho*u + sqrt(1-rho^2)*v
                sigma_pcn = min(lambda_pcn*sigma0,1);
                rho   = sqrt(1-sigma_pcn.^2);
                n_acp = 0;
                dN_cal = 0;
                idx_sam = idx_grp(iter,1):idx_grp(iter,2);

                for k=1:Ns
                    if k == 1 && sedsav
                        % start from seeds
                        u(:,idx_sam,k) = u_seed(:,idx_sam);
                        y(:,idx_sam,k) = y_seed(:,idx_sam);
                    else
                        % pCN proposal (Eqn. 22)
                        u(:,idx_sam,k) = rho.*u_seed(:,idx_sam) + sigma_pcn.*v(:,idx_sam,k);

                        % evaluate target
                        baystc = baystc.EvlY(u(:,idx_sam,k));
                        y(:,idx_sam,k) = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];

                        % MH accept/reject (Eqn. 23, log-scale)
                        idx_rej = exp(y(1,idx_sam,k)-y_seed(1,idx_sam)) < u_mh(1,idx_sam,k);
                        idx_sam_rej = idx_sam(idx_rej);

                        % revert rejected samples
                        u(:,idx_sam_rej,k) = u_seed(:,idx_sam_rej);
                        y(:,idx_sam_rej,k) = y_seed(:,idx_sam_rej);

                        % statistics
                        n_acp  = n_acp + length(idx_sam)-length(idx_sam_rej);
                        dN_cal = dN_cal + baystc.Ncal;

                        % update seeds
                        u_seed(:,idx_sam) = u(:,idx_sam,k);
                        y_seed(:,idx_sam) = y(:,idx_sam,k);
                    end
                end

                % acceptance rate
                acp = n_acp/length(idx_sam)/Ns;

                % adaptive scaling update
                lambda_pcn = exp(log(lambda_pcn)+(acp-acp_opt)/sqrt(iter));

                % update total calls
                N_cal = N_cal+dN_cal;
            end
        end
    end
end
