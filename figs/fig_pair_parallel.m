% fig_ellipse.m — separation and angular velocity vs time for ellipse pairs,
% spanning both oblate (a/b < 1) and prolate (a/b > 1) aspect ratios.
% Companion to fig_pair_fields.m (same examples/pair/output data).

close all;

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'parallel');
P    = palette();

files = dir(fullfile(dataDir, 'a*_d*.dat'));
a_vals = zeros(1, numel(files));
b_vals = zeros(1, numel(files));
for k = 1:numel(files)
    tok = sscanf(files(k).name, 'a%f_b%f_d%f.dat');
    a_vals(k) = tok(1);
    b_vals(k) = tok(2);
end
[~, isort] = sort(a_vals);
files  = files(isort);
a_vals = a_vals(isort);
b_vals = b_vals(isort);

cmap = cmocean('tempo', numel(files) + 2);
cmap = cmap(2:end-1, :);
cmap = colorcet('D11', 'N', numel(files)).^1.5;

fig   = journal_figure([6.75, 2.25], 2, 'tex');
axSep = axes(fig); hold(axSep, 'on'); box(axSep, 'on');
axOmg = axes(fig); hold(axOmg, 'on'); box(axOmg, 'on');

for k = 1:numel(files)
    A    = readmatrix(fullfile(files(k).folder, files(k).name));
    A    = A(all(~isnan(A), 2), :);
    c    = cmap(k, :);
    lbl  = ['{\it a/b} = ' num2str(a_vals(k)/b_vals(k), '%.2f')];
    t    = A(:,1);
    sep  = A(:,8);
    ang1 = unwrap(A(:,4));
    ang2 = unwrap(A(:,7));
    plot(axSep, t, sep/sep(1), 'Color', c, 'DisplayName', lbl);
    plot(axOmg, t, ang1 - pi/2, 'Color', c, 'HandleVisibility', 'off');
    plot(axOmg, t, ang2 - pi/2, 'Color', c, 'LineStyle', '--', 'HandleVisibility', 'off');

    if a_vals(k) > b_vals(k)
        plot(axSep, t(end), sep(end)/sep(1), 'o', 'Color', c, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
        plot(axOmg, t(end), ang1(end) - pi/2, 'o', 'Color', c, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
        plot(axOmg, t(end), ang2(end) - pi/2, 'o', 'Color', c, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    end
end

linkaxes([axSep, axOmg], 'x');
xlabel(axSep, 'time, {\it t}');
xlabel(axOmg, 'time, {\it t}');
ylabel(axSep, 'dist({\it\Gamma}_1, {\it\Gamma}_2)');
ylabel(axOmg, '{\it \theta} - \pi/2');
legend(axSep, 'Location', 'southwest', 'Color', 'none');

xlim(axSep, [0 500]); xlim(axOmg, [0 500]);
% ── Layout ────────────────────────────────────────────────────────────────────
ml = 1.0; mr = 0.20; mb = 0.75; mt = 0.30; gapH = 1.0;

fig.Units = 'inches'; figW = fig.Position(3); figH = fig.Position(4);
set([axSep, axOmg], 'Units', 'inches');
panelW = 0.975 * (figW - ml - mr - gapH) / 2;
panelH = figH - mb - mt;

axSep.Position = [ml,                mb, panelW, panelH];
axOmg.Position = [ml + panelW + gapH, mb, panelW, panelH];

panel_label('\it (a)', axSep, -0.95, 0.2);
panel_label('\it (b)', axOmg, -0.95, 0.2);

function a = panel_label(str, ax, xOffset, yOffset)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', ax.FontSize, 'Units', 'inches');
    a.Position(1) = ax.Position(1) + xOffset;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + yOffset;
end
