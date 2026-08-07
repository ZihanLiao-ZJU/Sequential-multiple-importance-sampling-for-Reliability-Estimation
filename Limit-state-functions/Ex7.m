classdef Ex7
    % Example 7: noisy PF (6D)   [paper Table 1]
    %
    % Definition
    % ---------------------------------------------------------------------
    % g(u) = u1 + 2*u2 + 2*u3 + u4 - 5*u5 - 5*u6 + 0.001 * sum_i sin(100*u_i)
    % Failure event:  g(u) < 0
    %
    % Configuration
    % ---------------------------------------------------------------------
    % 1) Dimension: 6
    % 2) Input RVs: independent lognormal distributions (mean/COV specification)
    % 3) Reference failure probability: Pf_ref = 5.29e-4
    %
    % Reference
    % ---------------------------------------------------------------------
    % [44] S. Engelund, R. Rackwitz, A benchmark study on importance sampling techniques in structural reliability, Structural Safety 12 (1993) 255-276. https://doi.org/10.1016/0167-4730(93)90056-7.
    %
    % Syntax
    % ---------------------------------------------------------------------
    % obj = Ex7
    properties
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Nfun       length of function list                          double            [1,1]
        % Dist       Joint distribution (lognormal marginals + Nataf)  Nataf             [1,1]
        % Ndim       number of dimensions                             double            [1,1]
        % -----------------------------------------------------------------------------------
        Nfun = 1;
        Dist
        Ndim = 6;
    end

    methods
        function obj = Ex7
            dist = [];
            sigma = zeros(6,1);
            mu = zeros(6,1);
            m = [120; 120; 120; 120; 50; 40];
            cov  = [8/120; 8/120; 8/120; 8/120; 10/50; 8/40];
            for i = 1:6
                sigma(i,:) = log(1+cov(i)^2)^0.5;
                mu(i,:) = log(m(i)/exp(sigma(i)^2/2));
                dist = cat(1,dist,makedist("Lognormal","mu",mu(i),"sigma",sigma(i)));
            end
            obj.Dist = dist;
            obj.Dist = Nataf(dist,eye(6));
        end
        function L = EvlLSF(~,theta)
            % Evaluate the performance function at samples theta
            % ----------------------------------------------------------------------
            % SYNTAX:
            % L = EvlLSF(obj,theta)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : constructed class
            % theta : samples in standard normal space                        [Ndim,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % L     : PF values; failure corresponds to L < 0                [1,Nsam]
            % ----------------------------------------------------------------------
            % REFERENCES:
            % ----------------------------------------------------------------------
            % [1] Engelund, S., & Rackwitz, R. (1993).
            % ----------------------------------------------------------------------

            % initialization
            L =  theta(1,:)+2*theta(2,:)+2*theta(3,:)+theta(4,:)-5*theta(5,:)-5*theta(6,:)+0.001*sum(sin(100*theta),1);
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
            text(theta_fpt(1,1)+0.5,theta_fpt(2,1),['$[',sprintf('%5.2f',theta_fpt(1,1)),',',sprintf('%5.2f',theta_fpt(2,1)),']$'],'Interpreter','latex','Color',"magenta");
            % figure settings
            colorbar
            grid on
            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            xlabel('$\theta_1$','Interpreter','latex')
            ylabel('$\theta_2$','Interpreter','latex')
            title('$\it{Piecewise\ Linear}$','Interpreter','latex')

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
            title('$\it{Piecewise\ Linear}$','Interpreter','latex')
            view(155,41)
        end

    end
end
