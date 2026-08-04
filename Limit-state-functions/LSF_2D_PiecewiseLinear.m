classdef LSF_2D_PiecewiseLinear
    % Piecewise-linear limit state function (LSF) in standard normal space
    %
    % Definition
    % ---------------------------------------------------------------------
    % L(theta) = min( L1(theta1), L2(theta2) )
    % Failure event:  L(theta) < 0
    %
    % Each component Li is a 1D piecewise-linear function with a breakpoint:
    %   L1(theta1) = { a1 - a2*theta1,  theta1 > t1
    %               { a3 - a4*theta1,  theta1 <= t1
    %   L2(theta2) = { b1 - b2*theta2,  theta2 > t2
    %               { b3 - b4*theta2,  theta2 <= t2
    %
    % Here, the breakpoints are fixed in code as:
    %   t1 = 3.5  (for theta1) ,  t2 = 2.0  (for theta2)
    %
    % Configuration
    % ---------------------------------------------------------------------
    % 1) Dimension: 2
    % 2) Reference distribution: standard normal (independent), implemented
    %    via a Nataf object with identity correlation
    % 3) Reported reference failure probability (for this parameter set):
    %    Pf_ref = 3.20e-5
    %
    % Reference
    % ---------------------------------------------------------------------
    % [1] Breitung, K. (2019). The geometry of limit state function graphs
    %     and subset simulation: Counterexamples. Reliability Engineering &
    %     System Safety, 182, 98–106.

    properties
        % Parameters and model settings
        % -----------------------------------------------------------------
        % Name   Description                                                     Type/Size
        % -----------------------------------------------------------------
        % Parm   LSF parameters packed as [a1;a2;a3;a4;b1;b2;b3;b4]              double [8,1]
        %        - (a1,a2,a3,a4): two linear pieces for L1(theta1)
        %        - (b1,b2,b3,b4): two linear pieces for L2(theta2)
        % Nfun   Number of component LSFs (scalar LSF => 1)                      double [1,1]
        % Ndim   Input dimension                                                 double [1,1]
        % Dist   Input joint distribution in physical space (here: standard U)    Nataf  [1,1]
        % -----------------------------------------------------------------
        Parm = [4;1;0.85;0.1;0.5;0.1;2.3;1];
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
            % - The LSF is constructed as the minimum of two 1D piecewise-linear
            %   components, producing a non-smooth (kinked) failure boundary.
            % - Breakpoints are hard-coded as theta1=3.5 and theta2=2.0.

            % initialization
            parm = obj.Parm;
            Nsam = size(theta,2);
            a1 = parm(1);
            a2 = parm(2);
            a3 = parm(3);
            a4 = parm(4);
            b1 = parm(5);
            b2 = parm(6);
            b3 = parm(7);
            b4 = parm(8);

            % evaluation
            ind1 = theta(1,:)>3.5;
            L1 = zeros(1,Nsam);
            L2 = zeros(1,Nsam);
            L1(ind1) = a1-a2*theta(1,ind1);
            L1(~ind1) = a3-a4*theta(1,~ind1);
            ind2 = theta(2,:)>2;
            L2(ind2) = b1-b2*theta(2,ind2);
            L2(~ind2) = b3-b4*theta(2,~ind2);
            L = min(L1,L2);
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

            % evaluation on grid
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

            % --- Figure 1: contour view ---
            figure
            contourf(xx,yy,zz,[zmin:abs(zmin):0,zmax/10:zmax/10:zmax]);
            hold on
            % failure boundary (L=0)
            theta_fpt = contour(xx,yy,zz,[0,0],"LineColor","red","LineWidth",1.5);

            % pick one contour vertex with minimal squared distance (visual cue)
            fpt_dis = sum(theta_fpt.^2);
            [~,ind_fpt] = min(fpt_dis,[],"all");
            theta_fpt = theta_fpt(:,ind_fpt);

            scatter(theta_fpt(1,:),theta_fpt(2,:),"Marker","o","MarkerEdgeColor","magenta","LineWidth",3);
            hold on
            text(theta_fpt(1,1)+0.5,theta_fpt(2,1), ...
                ['$[',sprintf('%5.2f',theta_fpt(1,1)),',',sprintf('%5.2f',theta_fpt(2,1)),']$'], ...
                'Interpreter','latex','Color',"magenta");

            colorbar
            grid on
            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            xlabel('$\theta_1$','Interpreter','latex')
            ylabel('$\theta_2$','Interpreter','latex')
            title('$\it{Piecewise\ Linear}$','Interpreter','latex')

            % --- Figure 2: surface view ---
            figure
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

            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            xlabel('$\theta_1$','Interpreter','latex')
            ylabel('$\theta_2$','Interpreter','latex')
            zlabel('$LSF$','Interpreter','latex')
            title('$\it{Piecewise\ Linear}$','Interpreter','latex')
            view(155,41)
        end
    end
end