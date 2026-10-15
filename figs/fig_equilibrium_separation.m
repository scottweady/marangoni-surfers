% fig_equilibrium_separation.m — far-field prediction for the equilibrium
% separation of an anti-parallel pair, r = eta/(4 pi U)(1 - beta/2),
% as a function of aspect ratio a/b at fixed area (ab = 1).
% The pair is offset perpendicular to the propulsion axis, so the particles
% touch at r = 2b; the prediction is unphysical below that curve.
% Uses the exact single-particle constants from single_particle.m.

close all;

P = palette();

ratio = logspace(log10(1/128), log10(128), 25)';
[U, eta, beta] = single_particle(ratio);
r  = eta ./ (4*pi*U) .* (1 - beta/2);
rc = 2 ./ sqrt(ratio);

fig = journal_figure([6.75/2, 2.5], 2, 'tex');
ax  = axes(fig);

fill(ax, [ratio; flipud(ratio)], [rc; zeros(size(rc))], 0.9 * [1 1 1], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
hold on
plot(ax, ratio, r, '-o', 'Color', P.c1, 'MarkerSize', 4, 'MarkerFaceColor', 'w', ...
    'DisplayName', '{\it{d}}^*');
plot(ax, ratio, rc, ':', 'Color', 'k', 'DisplayName', '2{\it{b}}');
ax.XScale = 'log';
% ax.XTick = [0.125 0.25 0.5 1 2 4 8];
% ax.XTickLabel = {'1/8', '1/4', '1/2', '1', '2', '4', '8'};
% xlim(ax, [1/16 16]);
ylim(ax, [0 6]);

text(ax, 0.9, 0.9, 'feasible ({\it{d}}^* > 2{\it{b}})', ...
     'Units', 'normalized', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'Margin', 2, ...
     'FontSize', ax.FontSize);

xlabel(ax, 'aspect ratio, {\it a/b}');
ylabel(ax, 'equilibrium separation, {\it{d}}^*');
% legend(ax, 'Location', 'northeast', 'Color', 'none');

% ── Layout ────────────────────────────────────────────────────────────────────
ml = 1.0; mr = 0.25; mb = 0.9; mt = 0.25;

fig.Units = 'inches'; figW = fig.Position(3); figH = fig.Position(4);
ax.Units = 'inches';
ax.Position = [ml, mb, figW - ml - mr, figH - mb - mt];
