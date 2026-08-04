classdef LSF_2D_LinearLog
    % Linear-log LSF in standard normal space
    %
    % Definition
    % ---------------------------------------------------------------------
    % L(theta) = min( L1(theta1), L2(theta2) )
    % Failure event:  L(theta) < 0
    %
    % L1(theta1) = a1 - a2*theta1
    % L2(theta2) = b1/(b2 + exp(-b3*(b4*theta2 + b5))) - c
    %
    % Configuration
    % ---------------------------------------------------------------------
    % 1. Dimension: 2
    % 2. Pf_ref = 3.23e-5
    % ---------------------------------------------------------------------
    % Reference:
    % ---------------------------------------------------------------------
    % [1] Breitung, K. The geometry of limit state function graphs and subset simulation:
    %     Counterexamples. Reliability Engineering & System Safety, 182 (2019), 98–106.
    %     -- Section 5.3
    % ---------------------------------------------------------------------
    properties
        % Parameters and model settings
        % -----------------------------------------------------------------
        % Name   Description                                                     Type/Size
        % -----------------------------------------------------------------
        % Parm   LSF parameters packed as [a1;a2;b1;b2;b3;b4;b5;c]               double [8,1]
        %        - (a1,a2): linear component for L1(theta1)
        %        - (b1..b5,c): logistic-type component for L2(theta2)
        % Nfun   Number of component LSFs (scalar LSF => 1)                      double [1,1]
        % Ndim   Input dimension                                                 double [1,1]
        % Dist   Input joint distribution in physical space (here: standard U)    Nataf [1,1]
        % -----------------------------------------------------------------
        Parm = [5;1;1;1;2;1;4;0.5];
        Nfun = 1;
        Ndim = 2;
        Dist = Nataf(repmat(makedist("Normal","mu",0,"sigma",1),2,1),eye(2));
    end

    methods
        function L = EvlLSF(obj,theta)
            % Evaluate the limit state function at samples theta
            %
            % Syntax
            % -----------------------------------------------------------------
            % L = obj.EvlLSF(theta)
            %
            % Inputs
            % -----------------------------------------------------------------
            % theta : samples in standard normal space                           [Ndim, Nsam]
            %
            % Outputs
            % -----------------------------------------------------------------
            % L     : LSF values; failure corresponds to L < 0                   [1,   Nsam]
            %
            % Notes
            % -----------------------------------------------------------------
            % - The LSF is the minimum of a linear term in theta1 and a smooth
            %   logistic-type term in theta2, leading to a non-smooth failure
            %   boundary where L1 = L2.

            % initialization
            parm = obj.Parm;
            a1 = parm(1);
            a2 = parm(2);
            b1 = parm(3);
            b2 = parm(4);
            b3 = parm(5);
            b4 = parm(6);
            b5 = parm(7);
            c  = parm(8);
           

            % evaluation
            L1 = a1 - a2*theta(1,:);
            L2 = b1./(b2+exp(-b3*(b4*theta(2,:)+b5)))-c;
            L  = min(L1,L2);
        end

        function PltFun(obj)
            % Visualize the LSF in 2D (contour + surface)
            %
            % Figures
            % -----------------------------------------------------------------
            % (1) Filled contour plot of L(theta) with the zero-level set L=0
            %     highlighted (failure boundary).
            % (2) Surface plot showing safe region (L>=0) and failure region
            %     (L<0), together with the plane L=0.
            %
            % Notes / limitations
            % -----------------------------------------------------------------
            % - The point marked as "failure point" is obtained from the contour
            %   output of MATLAB and is intended as a visual indicator (not a
            %   robust FORM/HL-RF design point solver).
            % - Grid evaluation uses a fixed range [-8,8]^2 with step 0.02.

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
            text(theta_fpt(1,1)-3.2,theta_fpt(2,1)-0.5,['$[',sprintf('%5.2f',theta_fpt(1,1)),',',sprintf('%5.2f',theta_fpt(2,1)),']$'],'Interpreter','latex','Color',"magenta");
            % figure settings
            colorbar    
            grid on
            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            xlabel('$\theta_1$','Interpreter','latex')
            ylabel('$\theta_2$','Interpreter','latex')
            title('$\it{Linear\ Log}$','Interpreter','latex')

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
            title('$\it{Linear\ Log}$','Interpreter','latex')
            view(-13,27)
        end

    end
end