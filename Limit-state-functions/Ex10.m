classdef Ex10
    % Example 10: linear hyperplane PF (nD)   [paper Table 4]
    %
    % Definition
    % ---------------------------------------------------------------------
    % g(u) = beta*sqrt(n) - sum(u)   (equivalently g(u) = beta - alpha'*u, ||alpha|| = 1)
    % default beta = 3;  Pf = 1.35e-3   (also provided: beta = 2,4,5,6 in the code)
    % Failure event:  g(u) < 0
    %
    % Configuration
    % ---------------------------------------------------------------------
    % 1) Dimension: n (user-specified)
    % 2) Reference distribution: independent standard normal
    % 3) Reference failure probability (beta = 3): Pf_ref = 1.35e-3
    %
    % Reference
    % ---------------------------------------------------------------------
    % [44] S. Engelund, R. Rackwitz, A benchmark study on importance sampling techniques in structural reliability, Structural Safety 12 (1993) 255-276. https://doi.org/10.1016/0167-4730(93)90056-7.
    %
    % Syntax
    % ---------------------------------------------------------------------
    % obj = Ex10(Ndim)
    properties(Constant)
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Nfun       length of function list                          double            [1,1]
        % -----------------------------------------------------------------------------------
        Nfun = 1;
    end

    properties
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Dist       Reference distribution (standard normal space)    Nataf             [1,1]
        % Ndim       number of dimensions                             double            [1,1]
        % -----------------------------------------------------------------------------------
        Dist
        Ndim
    end

    methods
        function obj = Ex10(Ndim)
            obj.Ndim = Ndim;
            obj.Dist = Nataf(repmat(makedist("Normal","mu",0,"sigma",1),Ndim,1),eye(Ndim));
        end

        function L = EvlLSF(~,theta)
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
            % [1] Dai & Wang (2009); Engelund & Rackwitz (1993).
            % ----------------------------------------------------------------------

            % initialization
            beta = 3;
            ndim = size(theta,1);

            % evaluation
            L = beta*sqrt(ndim) - sum(theta,1);
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
            text(theta_fpt(1,1)+0.5,theta_fpt(2,1),['$[',sprintf('%5.2f',theta_fpt(1,1)),',',sprintf('%5.2f',theta_fpt(2,1)),']$'],'Interpreter','latex','Color',"magenta");
            % figure settings
            colorbar
            grid on
            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            xlabel('$\theta_1$','Interpreter','latex')
            ylabel('$\theta_2$','Interpreter','latex')
%             title('$\it{Linear\ PF}$','Interpreter','latex')

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
%             title('$\it{Linear\ PF}$','Interpreter','latex')
            view(67,22)
        end

    end
end
