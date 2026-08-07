classdef Ex2
    % Example 2: PF with changing topology (2D)   [paper Table 1]
    %
    % Definition
    % ---------------------------------------------------------------------
    % g(u) = d + c1 / ( ( b1*(u1-a1)^2 + b2*(u2-a2)^2 )^2 + 1 )
    %           + c2 / ( ( b3*(u1-a3)^2 + b4*(u2-a4)^2 )^2 + 1 )
    % with (a1,a2) = (-2,0), (a3,a4) = (2.5,0.5), (b1,b2) = (4/9,1/25),
    %      (b3,b4) = (1/4,1/25), c1 = 30, c2 = 20, d = -5
    % Failure event:  g(u) < 0
    %
    % Configuration
    % ---------------------------------------------------------------------
    % 1) Dimension: 2
    % 2) Reference distribution: independent standard normal
    % 3) Reference failure probability: Pf_ref = 1.13e-5
    %
    % Reference
    % ---------------------------------------------------------------------
    % [16] K. Breitung, The geometry of limit state function graphs and subset simulation: Counterexamples, Reliability Engineering & System Safety 182 (2019) 98-106. https://doi.org/10.1016/j.ress.2018.10.008.
    %
    % Syntax
    % ---------------------------------------------------------------------
    % obj = Ex2
    properties
        % Parameters and model settings
        % -----------------------------------------------------------------
        % Name   Description                                                     Type/Size
        % -----------------------------------------------------------------
        % Parm   PF parameters packed as [a1;a2;a3;a4;b1;b2;b3;b4;c1;c2;d]      double [11,1]
        %        - (a1,a2) and (a3,a4): centres in (theta1,theta2) for two terms
        %        - (b1,b2) and (b3,b4): anisotropic weights for the squared radii
        %        - (c1,c2): amplitudes; d: vertical shift                              
        % Nfun   Number of component LSFs (scalar PF => 1)                      double [1,1]
        % Ndim   Input dimension                                                 double [1,1]
        % Dist   Input joint distribution in physical space (here: standard U)    Nataf [1,1]
        % -----------------------------------------------------------------
        Parm = [-2;0;2.5;0.5;4/9;1/25;1/4;1/25;30;20;-5];
        Nfun = 1;
        Ndim = 2;
        Dist = Nataf(repmat(makedist("Normal","mu",0,"sigma",1),2,1),eye(2));
    end

    methods
        function L = EvlLSF(obj,theta)
            % Evaluate the performance function at samples theta
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
            % L     : PF values; failure corresponds to L < 0                   [1,   Nsam]
            %
            % Notes
            % -----------------------------------------------------------------
            % - The PF is a superposition of two smooth, localized terms whose
            %   interaction may produce a failure boundary with changing topology
            %   as parameters vary.

            % initialization
            parm = obj.Parm;
            a1 = parm(1);
            a2 = parm(2);
            a3 = parm(3);
            a4 = parm(4);
            b1 = parm(5);
            b2 = parm(6);
            b3 = parm(7);
            b4 = parm(8);
            c1 = parm(9);
            c2 = parm(10);
            d  = parm(11);

            % evaluation
            L = d + c1./((b1*(theta(1,:)-a1).^2+b2*(theta(2,:)-a2).^2).^2+1) + c2./((b3*(theta(1,:)-a3).^2+b4*(theta(2,:)-a4).^2).^2+1);
        end

        function PltFun(obj)
            % Visualize the PF in 2D (contour + surface)
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
            text(theta_fpt(1,1)-3.2,theta_fpt(2,1),['$[',sprintf('%5.2f',theta_fpt(1,1)),',',sprintf('%5.2f',theta_fpt(2,1)),']$'],'Interpreter','latex','Color',"magenta");
            % figure settings
            colorbar
            grid on
            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            xlabel('$\theta_1$','Interpreter','latex')
            ylabel('$\theta_2$','Interpreter','latex')
%             title('$\it{Changing\ Toplogical}$','Interpreter','latex')

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
%             title('$\it{Changing\ Toplogical}$','Interpreter','latex')
            % title('$\it{Changing\ Toplogical}$','Interpreter','latex')
            view(-13,27)
        end

    end
end
