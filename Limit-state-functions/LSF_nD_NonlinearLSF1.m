classdef LSF_nD_NonlinearLSF1
    % Nonlinear LSF in standard normal space
    % Syntax
    % ---------------------------------------------------------------------
    % obj = LSF_nD_NonlinearLSF1(Ndim)
    % ---------------------------------------------------------------------
    % Configuration
    % ---------------------------------------------------------------------
    % 1. Dimension: n (user-specified)
    % 2. Pf_ref = 4.75e-6; (k = -10.0, default)
    %    Pf = 6.41e-5;     (k = 0.2)
    %    Pf = 1.40e-3;     (k = 0.6)
    %    Pf = 8.90e-3;     (k = 1.0)
    %    Pf = 1.39e-5;     (k = -1.0)
    %    Pf = 6.50e-6;     (k = -5.0)
    % ---------------------------------------------------------------------
    % Reference:
    % ---------------------------------------------------------------------
    % [1] Ghalehnovi, M., Rashki, M., & Ameryan, A. First order control variates
    %     algorithm for reliability analysis of engineering structures. Applied
    %     Mathematical Modelling, 77 (2020), 829–847.
    %     -- Section 5.2, Case 2
    % ---------------------------------------------------------------------


    properties(Constant)
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Nfun       length of function list                          double            [1,1]
        % Pf_ref     reference failure probability (default setting)  double            [1,1]
        % -----------------------------------------------------------------------------------
        Nfun = 1;
    end

    properties
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Ndim       number of dimensions                             double            [1,1]
        % Dist       Reference distribution (standard normal space)    Nataf             [1,1]
        % -----------------------------------------------------------------------------------
        Ndim
        Dist
    end

    methods
        function obj = LSF_nD_NonlinearLSF1(Ndim)
            obj.Ndim = Ndim;
            obj.Dist = Nataf(repmat(makedist("Normal","mu",0,"sigma",1),Ndim,1),eye(Ndim));
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
            % L     : function value; failure corresponds to L < 0                [1,Nsam]
            % ----------------------------------------------------------------------
            % REFERENCES:
            % ----------------------------------------------------------------------
            % [1] Ghalehnovi et al. (2020).
            % ----------------------------------------------------------------------

            % initialization
            beta0 = 4;
            k = -10;

            % evaluation
            L = beta0 - sum(theta,1)/sqrt(obj.Ndim) -  k*(theta(1,:)-theta(2,:)).^2/4;
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
            % evlauation
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
%             title('$\it{Nonlinear\ LSF}$','Interpreter','latex')

            figure
            % plot the LSF surf
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
            zlabel('$LSF$','Interpreter','latex')
%             title('$\it{Nonlinear\ LSF}$','Interpreter','latex')
            view(155,41)
        end

    end
end