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
        function [u,y,Ncal,baystc] = SamGen(obj,u_sed,y_sed,Ns)
            % Generate MCMC samples targeting prior × likelihood
            % pCN proposal + adaptive scaling on acceptance rate

            % initialization
            baystc = obj.BayStc;
            sedsav = obj.SedSav;
            Ndim = obj.BayStc.Ndim;
            Nfun = obj.BayStc.Nfun;
            Nite = obj.NumIte;
            lambda = obj.Lambda;
            acpopt = obj.AcpOpt;

            % linear map for correlated Gaussian proposal
            if isfield(baystc.G,'L')
                L = baystc.G.L;
            else
                L = eye(Ndim);
            end

            Nc = size(u_sed,2);
            Ncal = 0;

            % allocate memory
            u = zeros(Ndim,Nc,Ns);
            y = zeros(Nfun,Nc,Ns);

            % Gaussian innovations and MH uniforms
            du = pagemtimes(L,randn(Ndim,Nc,Ns));
            p  = rand(1,Nc,Ns);

            % group chains for adaptive updates
            Ngrp = max(ceil(Nc/Nite),1);
            Nite = ceil(Nc/Ngrp);
            ind_grp = min(repmat([1,Ngrp],Nite,1)+(0:Ngrp:(Nite-1)*Ngrp)',Nc);

            % base proposal scale
            % sigma0 = max(std(u_sed,0,2),eps);
            sigma0 = 1;

            for iter = 1:Nite
                % pCN parameters
                sigma = min(lambda*sigma0,1);
                rho   = sqrt(1-sigma.^2);
                Nacp  = 0;
                dNcal = 0;
                ind_sam = ind_grp(iter,1):ind_grp(iter,2);

                for k=1:Ns
                    if k == 1 && sedsav
                        % start from seeds
                        u(:,ind_sam,k) = u_sed(:,ind_sam);
                        y(:,ind_sam,k) = y_sed(:,ind_sam);
                    else
                        % pCN proposal
                        u(:,ind_sam,k) = rho.*u_sed(:,ind_sam) + sigma.*du(:,ind_sam,k);

                        % evaluate target
                        baystc = baystc.EvlY(u(:,ind_sam,k));
                        y(:,ind_sam,k) = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];

                        % MH accept/reject (log-scale)
                        ind_rej = exp(y(1,ind_sam,k)-y_sed(1,ind_sam)) < p(1,ind_sam,k);
                        ind_sam_rej = ind_sam(ind_rej);

                        % revert rejected samples
                        u(:,ind_sam_rej,k) = u_sed(:,ind_sam_rej);
                        y(:,ind_sam_rej,k) = y_sed(:,ind_sam_rej);

                        % statistics
                        Nacp  = Nacp + length(ind_sam)-length(ind_sam_rej);
                        dNcal = dNcal + baystc.Ncal;

                        % update seeds
                        u_sed(:,ind_sam) = u(:,ind_sam,k);
                        y_sed(:,ind_sam) = y(:,ind_sam,k);
                    end
                end

                % acceptance rate
                acp = Nacp/length(ind_sam)/Ns;

                % adaptive scaling update
                lambda = exp(log(lambda)+(acp-acpopt)/sqrt(iter));

                % update total calls
                Ncal = Ncal+dNcal;
            end
        end
    end
end
