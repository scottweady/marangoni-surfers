% fig_stability.m — maximal growth rate of the linearized mean-field model
% about the von Mises state psi_0 ∝ exp(kappa cos theta), as a function of
% wavenumber k and wavevector angle phi, for kappa = 5, 10, 15.
% Insets show the base distribution psi_0(theta).
%
% Requires: data/stability/growth_rate.mat (run scripts/meanfield/stability_solve.m)

close all

load(fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'stability', 'growth_rate.mat'));
[K, PH] = meshgrid(kh, ph);
smax = max(cellfun(@(S) max(S(:)), Smax));

P = palette();
cmap = P.cmap();
cmap = [1 1 1; cmap];

pc = linspace(1e-3, pi/2, 200);
kc = beta*(1 - lambda*sin(pc).^2).^2 ./ (8*lambda*U*sin(pc).^2);

t = linspace(-pi, pi, 400);
psi_0 = @(kappa) exp(kappa*(cos(t) - 1)) / (2*pi*besseli(0, kappa, 1));

fig = journal_figure([6.75 2.5], 2, 'tex');
labels = {'\it (a)', '\it (b)', '\it (c)'};
for j = 1:3
    ax(j) = axes(fig, 'Layer', 'top'); hold(ax(j), 'on'); box(ax(j), 'on');
    S = Smax{j};
    S(S<1e-3) = NaN;
    keep = any(~isnan(S), 2);   % skip rows with no unstable points
    contourf(ax(j), K(keep,:), PH(keep,:), S(keep,:), linspace(0, smax, 9), 'LineColor', 'none', 'HandleVisibility', 'off');
    plot(ax(j), kc, pc, 'k:', 'LineWidth', 2, 'DisplayName', '{\it k_c} ({\it \kappa} \rightarrow \infty)');
    colormap(ax(j), cmap); clim(ax(j), [0 smax]);
    xlim(ax(j), kh([1 end])); ylim(ax(j), ph([1 end]));
    yticks(ax(j), [0 pi/8 pi/4 3*pi/8 pi/2]);
    yticklabels(ax(j), {'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'});
    xlabel(ax(j), '{\it{k}^\prime}');
end
ylabel(ax(1), '{\it{\phi}}');
yticklabels(ax(2:3), {});
legend(ax(1), 'Location', 'southwest', 'Color', 'none');

cb = colorbar(ax(3));
cb.Title.String = 'max Re({\it\sigma}^\prime)';
cb.Title.FontName = 'Times';
cb.Title.FontSize = ax(3).FontSize;

% ── Layout ────────────────────────────────────────────────────────────────────
ml = 1.1; mr = 1.3; mb = 0.7; mt = 0.4; gap = 0.7;

fig.Units = 'inches'; figW = fig.Position(3); figH = fig.Position(4);
w = (figW - ml - mr - 2*gap)/3; h = figH - mb - mt;
for j = 1:3
    ax(j).Units = 'inches';
    ax(j).Position = [ml + (j-1)*(w + gap), mb, w, h];
    panel_label(labels{j}, ax(j), -0.95 + 0.3*(j > 1), 0.1);
end
cb.Units = 'inches';
cb.Position = [ml + 3*w + 2*gap + 0.12, mb, 0.12, h];

% ── Insets: base distribution ─────────────────────────────────────────────────
wi = 0.42*w; hi = 0.32*h;
for j = 1:3
    axi = axes(fig, 'Units', 'inches', 'Position', [sum(ax(j).Position([1 3])) - wi - 0.1, sum(ax(j).Position([2 4])) - hi - 0.1, wi, hi]);
    hold(axi, 'on'); box(axi, 'on');
    plot(axi, t, psi_0(kappas(j)), 'k');
    xlim(axi, [-pi pi]); ylim(axi, [0 1.1*max(psi_0(max(kappas)))]);
    xticks(axi, [-pi 0 pi]); xticklabels(axi, {'-\pi', '0', '\pi'});
    yticks(axi, []);
    axi.FontSize = 0.7*ax(j).FontSize;
    xlabel(axi, '{\it \theta}');
    ylabel(axi, '{\it \psi_0}');
    text(axi, 0.95, 0.9, sprintf('{\\it \\kappa} = %g', kappas(j)), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', axi.FontSize);
end

function a = panel_label(str, ax, xOffset, yOffset)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', 20, 'Units', 'inches');
    a.Position(1) = ax.Position(1) + xOffset;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + yOffset;
end
