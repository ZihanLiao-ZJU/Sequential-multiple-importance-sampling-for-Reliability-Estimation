clear; clc; close all;

%% problem
Ndim = 2;
% Tfun = LSF_nD_LinearLSF1(Ndim);
Tfun = LSF_nD_NonlinearLSF1(Ndim);

%% SeMIS-SuS
IntBay = SeMIS_SuS(Tfun);
sampler = aCS;

semcs = SeMIS(IntBay, sampler);
semcs.NumSam = 10000;

out_semcs = semcs.RunIte;

%% plot
plot_SeMIS_SuS_iterations(out_semcs, IntBay);


function plot_SeMIS_SuS_iterations(out, IntBay)
    Nite = out.Nite;

    all_u = [];
    for ite = 1:Nite
        if ~isempty(out.x{ite})
            all_u = [all_u, out.x{ite}(1:end,:)];
        end
    end

    if isempty(all_u)
        error('No samples found.');
    end

    u1_min = min(all_u(1,:));
    u1_max = max(all_u(1,:));
    u2_min = min(all_u(2,:));
    u2_max = max(all_u(2,:));

    du1 = max(u1_max - u1_min, 1);
    du2 = max(u2_max - u2_min, 1);

    xlim_all = [u1_min - 0.15*du1, u1_max + 0.15*du1];
    ylim_all = [u2_min - 0.15*du2, u2_max + 0.15*du2];

    ncol = min(3, Nite);
    nrow = ceil(Nite / ncol);

    figure('Color','w','Position',[100,100,420*ncol,360*nrow]);
    tiledlayout(nrow, ncol, 'TileSpacing','compact', 'Padding','compact');

    for ite = 1:Nite
        nexttile; hold on; box on; grid on;

        u = out.x{ite}(1:end,:);

        h_line_all = gobjects(0);

        % -------------------- Samples --------------------
        h_samples = scatter(u(1,:), u(2,:), 18, ...
            'MarkerFaceColor',[0.75 0.75 0.75], ...
            'MarkerEdgeColor','none');

        % -------------------- Hyperplanes --------------------
        if ite + 1 <= numel(out.g) && ~isempty(out.g{ite+1})
            g_next = out.g{ite+1};

            if isfield(g_next,'Beta') && ~isempty(g_next.Beta) && ...
               isfield(g_next,'f_hat_th') && ~isempty(g_next.f_hat_th)

                Nhp = size(g_next.Beta, 2);
                cmap = lines(Nhp);

                for m = 1:Nhp
                    Beta_m = g_next.Beta(:,m);
                    lhat_m = g_next.f_hat_th(m);

                    h_line = plot_pls_line(Beta_m, lhat_m, xlim_all, ylim_all);
                    set(h_line, ...
                        'Color', cmap(m,:), ...
                        'LineWidth', 1.8);

                    h_line_all(end+1) = h_line; %#ok<AGROW>
                end

                % -------------------- Seeds --------------------
                intBay_tmp = IntBay;
                intBay_tmp.G = g_next;

                if ite < Nite
                    Nseed_target = out.Nsam(ite+1);
                else
                    Nseed_target = size(u,2);
                end

                [x_sed, ~] = SedSlt(intBay_tmp, out.x{ite}, out.y{ite}, Nseed_target);
                x_sed = x_sed(1:end,:);

                if ~isempty(x_sed)
                    h_seeds = scatter(x_sed(1,:), x_sed(2,:), 28, ...
                        'o', ...
                        'MarkerEdgeColor', [0.20 0.20 0.90], ...
                        'MarkerFaceColor', 'none', ...
                        'LineWidth', 1.1);
                else
                    h_seeds = gobjects(1);
                end

            else
                h_seeds = gobjects(1);
            end
        else
            h_seeds = gobjects(1);
        end

        % -------------------- Labels --------------------
        xlabel('$u_1$', ...
            'Interpreter', 'latex', ...
            'FontName', 'Times New Roman', ...
            'FontSize', 16);

        ylabel('$u_2$', ...
            'Interpreter', 'latex', ...
            'FontName', 'Times New Roman', ...
            'FontSize', 16);

        title(['Iteration $i = ', num2str(ite), '$'], ...
            'Interpreter', 'latex', ...
            'FontName', 'Times New Roman', ...
            'FontSize', 16);

        set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

        xlim(xlim_all);
        ylim(ylim_all);
        axis equal;

        % -------------------- Legend: Samples -> Seeds -> Hyperplanes --------------------
        lgd_handles = h_samples;
        lgd_names   = {'Samples'};

        if exist('h_seeds','var') && isgraphics(h_seeds)
            lgd_handles(end+1) = h_seeds;
            lgd_names{end+1} = 'Seeds';
        end

        for m = 1:numel(h_line_all)
            lgd_handles(end+1) = h_line_all(m);
            lgd_names{end+1} = sprintf('Hyperplane $m=%d$', m);
        end

        legend(lgd_handles, lgd_names, ...
            'Location', 'best', ...
            'FontName', 'Times New Roman', ...
            'FontSize', 16, ...
            'Interpreter', 'latex');
    end
end


function h = plot_pls_line(Beta, fhat_th, xlim_all, ylim_all)
    b0 = Beta(1);
    b1 = Beta(2);
    b2 = Beta(3);

    if abs(b2) > 1e-12
        u1 = linspace(xlim_all(1), xlim_all(2), 400);
        u2 = (fhat_th - b0 - b1*u1) / b2;
        ind = (u2 >= ylim_all(1)) & (u2 <= ylim_all(2));

        if any(ind)
            h = plot(u1(ind), u2(ind), '-');
        else
            h = plot(u1, u2, '-');
        end

    elseif abs(b1) > 1e-12
        u1c = (fhat_th - b0) / b1;
        h = plot([u1c, u1c], ylim_all, '-');

    else
        h = plot(nan, nan, '-');
    end
end