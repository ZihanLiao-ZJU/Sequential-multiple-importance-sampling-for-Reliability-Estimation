classdef Ex9
    % Example 9: two-degree-of-freedom damped oscillator PF (8D)   [paper Table 1]
    %
    % Definition
    % ---------------------------------------------------------------------
    % g(u) = Fs - a1*ks*sqrt( pi*S0/(a2*zetas*omegas^3) * ( ... ) )
    % (displacement response of the 2-DOF damped oscillator; full expression in Ref. [46])
    % Failure event:  g(u) < 0
    %
    % Configuration
    % ---------------------------------------------------------------------
    % 1) Dimension: 8
    % 2) Input RVs: independent lognormal distributions (mean/COV specification)
    % 3) Reference failure probability: Pf_ref = 3.78e-7
    %
    % Reference
    % ---------------------------------------------------------------------
    % [46] A. Der Kiureghian, M. De Stefano, Efficient algorithm for second-order reliability analysis, Journal of Engineering Mechanics 117 (1991) 2904-2923. https://doi.org/10.1061/(ASCE)0733-9399(1991)117:12(2904).
    %
    % Syntax
    % ---------------------------------------------------------------------
    % obj = Ex9
    properties(Constant)
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Parm       parameters                                       double         [Npar,1]
        % Nfun       length of function list                          double            [1,1]
        % Ndim       number of dimensions                             double            [1,1]
        % -----------------------------------------------------------------------------------
        Parm = [3;4;4;4];
        Nfun = 1;
        Ndim = 8;
    end

    properties
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Dist       Joint distribution (lognormal marginals + Nataf)  Nataf             [1,1]
        % -----------------------------------------------------------------------------------
        Dist
    end

    methods
        function obj = Ex9
            Ndim = obj.Ndim;
            dist = [];
            sigma = zeros(Ndim,1);
            mu = zeros(Ndim,1);
            m = [1.5; 0.01; 1; 0.01; 0.05; 0.02; 27.5; 100];
            cov  = [0.1; 0.1; 0.2; 0.2; 0.4; 0.5; 0.1; 0.1];
            for i = 1:Ndim
                sigma(i,:) = log(1+cov(i)^2)^0.5;
                mu(i,:) = log(m(i)/exp(sigma(i)^2/2));
                dist = cat(1,dist,makedist("Lognormal","mu",mu(i),"sigma",sigma(i)));
            end
            obj.Dist = Nataf(dist,eye(8));
        end

        function L = EvlLSF(obj,theta)
            % evaluate the function value
            % ----------------------------------------------------------------------
            % SYNTAX:
            % L = EvlLSF(obj,theta)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : constructed class
            % theta : variable samples                                   [Ndim,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % L     : function value; failure corresponds to L < 0                 [1,Nsam]
            % ----------------------------------------------------------------------
            % REFERENCES:
            % ----------------------------------------------------------------------
            % [1] Kiureghian, A. D., & Stefano, M. D. (1991).
            % ----------------------------------------------------------------------


            % evaluation
            parm = obj.Parm;
            a1 = parm(1);
            a2 = parm(2);
            a3 = parm(3);
            a4 = parm(4);
            omegas = sqrt(theta(4,:)./theta(2,:));
            omegap = sqrt(theta(3,:)./theta(1,:));
            omegaa = (omegas+omegap)/2;
            zetas = theta(6,:);
            zetap = theta(5,:);
            zetaa = (zetas+zetap)/2;
            Fs = theta(7,:);
            ks = theta(4,:);
            S0 = theta(8,:);
            gamma = theta(2,:)./theta(1,:);
            theta = (omegap-omegas)./omegaa;
            L = Fs - a1*ks.*sqrt(pi*S0/a2./zetas./omegas.^3.*(zetaa.*zetas.*(zetap.*omegap.^3+zetas.*omegas.^3).*omegap./(zetap.*zetas.*(a3*zetaa.^2+theta.^2)+gamma.*zetaa.^2)./(a4*zetaa.*omegaa.^4)));
        end

        function PltFun(obj)
            % plot the function
            % ----------------------------------------------------------------------
            % SYNTAX:
            % PltFun (obj)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj  : constructed class
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % ----------------------------------------------------------------------
            % REFERENCES:
            % ----------------------------------------------------------------------
            % [1].
            % ----------------------------------------------------------------------

            x_ref = -8:0.02:8;
            y_ref = -8:0.02:8;
            xmin = min(x_ref);
            xmax = max(x_ref);
            Nx = length(x_ref);
            ymin = min(y_ref);
            ymax = max(y_ref);
            Ny = length(y_ref);
            [xx,yy] = meshgrid(x_ref,y_ref);
            zz = zeros(Nx,Ny);
            % evaluation
            for iy = 1:Ny
                zz(iy,:) = obj.EvlLSF([xx(iy,:);yy(iy,:)]);
            end
            zmin = min(zz,[],"all");
            zmax = max(zz,[],"all");
            ind0 = zz<0;
            zz_s = zz;
            zz_s(ind0) = nan;
            zz_f = zz;
            zz_f(~ind0) = nan;

            % figure plot
            figure
            % plot the likelihood cloud and contour
            contourf(xx,yy,zz,[zmin:abs(zmin):0,zmax/10:zmax/10:zmax]);
            hold on
            % the failure point
            theta_fpt = contour(xx,yy,zz,[0,0],"LineColor","red","LineWidth",1.5);
            fpt_dis = sum(theta_fpt.^2);
            [~,ind_fpt] = min(fpt_dis,[],"all");
            theta_fpt = theta_fpt(:,ind_fpt);
            % scatter the failure points
            scatter(theta_fpt(1,:),theta_fpt(2,:),"Marker","o","MarkerEdgeColor","magenta","LineWidth",3);
            hold on
            % label the mle points
            text(theta_fpt(1,1)-1.5,theta_fpt(2,1)-0.6,['$[',sprintf('%5.2f',theta_fpt(1,1)),',',sprintf('%5.2f',theta_fpt(2,1)),']$'],'Interpreter','latex','Color',"magenta");
            % figure settings
            colorbar
            grid on
            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            xlabel('$\theta_1$','Interpreter','latex')
            ylabel('$\theta_2$','Interpreter','latex')
            % title('$\it{Nonlinear\ PF}$','Interpreter','latex')

            figure
            % plot the PF surf
            surf(xx,yy,zz_s,"EdgeColor","none","FaceColor","green","FaceAlpha",0.3);
            hold on
            surf(xx,yy,zz_f,"EdgeColor","none","FaceColor","red");
            hold on
            surf(xx,yy,zeros(Nx,Ny),"EdgeColor","none","FaceColor","blue","FaceAlpha",0.3);
            hold on
            plot3(xx(1:20:Ny,:)',yy(1:20:Ny,:)',zz(1:20:Ny,:)',"Color","black","LineStyle","-","LineWidth",0.2)
            hold on
            plot3(xx(:,1:20:Nx),yy(:,1:20:Nx),zz(:,1:20:Nx),"Color","black","LineStyle","-","LineWidth",0.2)
            hold on
            scatter3(theta_fpt(1,:),theta_fpt(2,:),0,"Marker","o","MarkerEdgeColor","magenta","LineWidth",3);
            % settings
            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            xlabel('$\theta_1$','Interpreter','latex')
            ylabel('$\theta_2$','Interpreter','latex')
            zlabel('$PF$','Interpreter','latex')
            % title('$\it{Nonlinear\ PF}$','Interpreter','latex')
            view(155,41)
        end

    end
end
