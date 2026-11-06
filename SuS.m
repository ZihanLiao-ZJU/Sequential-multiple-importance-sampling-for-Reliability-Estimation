classdef SuS
    % Scalar subset-simulation ISD shared by the reliability toolboxes.
    % Both selectors sample/store Gaussian u. 'standard_normal' maps through
    % PF.Dist; the historical name 'physical' maps through PF.Rdis.

    properties
        TF
        Nfun
        Ndim
        Ncal = 0
        X
        Y
        G = struct('parameters', inf, 'p', 0.1)
        Nsam = 1000
        CoordinateSystem = 'standard_normal'
        ThresholdRule = 'hazen'
        ClipThreshold = true
    end

    properties (Dependent)
        p
        FlgCvg
    end

    methods
        function obj = SuS(PF, coordinate_system, threshold_rule, clip_threshold)
            obj.TF = PF;
            obj.Nfun = PF.Nfun + 3;
            obj.Ndim = PF.Ndim;
            obj.G.L = eye(obj.Ndim);
            if nargin >= 2
                obj.CoordinateSystem = coordinate_system;
            end
            if nargin >= 3
                obj.ThresholdRule = threshold_rule;
            end
            if nargin >= 4
                obj.ClipThreshold = clip_threshold;
            end
            assert(any(strcmp(obj.CoordinateSystem, {'standard_normal', 'physical'})));
            assert(any(strcmp(obj.ThresholdRule, {'hazen', 'midpoint'})));
        end

        function value = get.p(obj)
            value = obj.G.p;
        end

        function obj = set.p(obj, value)
            obj.G.p = value;
        end

        function obj = EvlY(obj, u)
            obj.X = u;
            if strcmp(obj.CoordinateSystem, 'physical')
                f = obj.TF.EvlLSF(obj.TF.Rdis.U2X(u));
                log_phi = -.5 * (obj.Ndim * log(2*pi) + sum(u.^2, 1));
            else
                f = obj.TF.EvlLSF(obj.TF.Dist.U2X(u));
                log_phi = logGauss(u);
            end
            obj.Y = [log_phi; f];
            obj.Ncal = size(u, 2);
        end

        function value = EvlLKF(obj)
            f = obj.Y(2, :);
            if strcmp(obj.CoordinateSystem, 'physical') && isinf(obj.G.parameters)
                value = zeros(size(f));
            else
                value = log(double(f < obj.G.parameters));
            end
        end

        function value = EvlPDF(obj, varargin)
            if ~isempty(varargin)
                value = -.5 * (log(2*pi) + obj.X(varargin{1}, :).^2);
            elseif ~isempty(obj.Y)
                value = obj.Y(1, :);
            else
                value = -.5 * (obj.Ndim * log(2*pi) + sum(obj.X.^2, 1));
            end
        end

        function [y, g, obj, N_cal] = UpdObj(obj, ~, y, p)
            if nargin < 4
                p = obj.G.p;
            end
            f = y(4, :);
            if strcmp(obj.ThresholdRule, 'midpoint')
                l = sort(f, 'ascend');
                n = round(numel(f) * p);
                g = (l(n) + l(n+1)) / 2;
            else
                g = quantile(f, p);
            end
            if obj.ClipThreshold
                g = max(g, 0);
            end
            obj.G.parameters = g;
            obj.G.p = p;
            g = obj.G;
            N_cal = 0;
        end

        function [x_seed, y_seed] = SltSed(obj, x, y, N)
            idx = y(4, :) < obj.G.parameters;
            n = sum(idx, 'all');
            if n <= 0
                x_seed = [];
                y_seed = [];
                return;
            end
            ns = round(N ./ (1:N));
            nc = round(N ./ ns);
            nc = nc(find(n-nc >= 0, 1, 'last'));
            sel = randperm(n, nc);
            xx = x(:, idx);
            yy = y(:, idx);
            x_seed = xx(:, sel);
            y_seed = yy(:, sel);
            if ~isempty(y_seed)
                obj.X = x_seed;
                obj.Y = y_seed(3:end, :);
                y_seed = [obj.EvlLKF; obj.EvlPDF; obj.Y];
            end
        end

        function x = SamGen(obj, N)
            x = randn(obj.Ndim, N);
        end

        function value = get.FlgCvg(obj)
            value = obj.G.parameters <= 0;
        end
    end
end
