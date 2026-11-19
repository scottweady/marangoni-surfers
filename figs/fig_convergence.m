% fig_convergence.m — convergence of the single-particle and pair solvers.
%
% Panel (a): single-particle error against the analytical center-of-mass velocity vs
%            boundary resolution N_t, for several aspect ratios a/b at fixed area.
% Panel (b): self-convergence of a close antiparallel pair vs expansion order M, boundary
%            held fixed; color is the aspect ratio, line style is the edge-to-edge gap.
%
% Requires: data/convergence/convergence_boundary.dat   (scripts/convergence/convergence.jl)
%           data/convergence/pair_gap_expansion.dat     (scripts/convergence/pair_gap_convergence.jl)

close all;

baseDir = fullfile(dataRoot, 'convergence');

[N, errN, hdrN] = load_errors(fullfile(baseDir, 'convergence_boundary.dat'));
[M, errM, hdrM] = load_errors(fullfile(baseDir, 'pair_gap_expansion.dat'));

ratios     = hdrN.ratios;
pairRatios = hdrM.ratios;
pairGaps   = hdrM.gaps;
gaps       = unique(pairGaps, 'stable');
styles     = {'-', '--', ':', '-.'};

% ── Colors ────────────────────────────────────────────────────────────────────
nr = numel(ratios);
P = palette();
cmap = P.linemap(nr + 2);
cmap = cmap(2:nr+1, :);
color_of = @(r) cmap(abs(ratios - r) < 1e-12, :);

% ── Figure ────────────────────────────────────────────────────────────────────
fig = journal_figure([6.75, 2.25], 2, 'tex');

axN = axes(fig); hold(axN, 'on'); box(axN, 'on');
axM = axes(fig); hold(axM, 'on'); box(axM, 'on');

for k = 1:nr
    plot(axN, N, errN(:,k), '-o', 'Color', cmap(k,:), 'MarkerSize', 4, ...
        'MarkerFaceColor', P.bg, 'DisplayName', ['{\it a/b} = ' ratio_label(ratios(k))]);
end

for j = 1:numel(pairRatios)
    ls = styles{find(gaps == pairGaps(j), 1)};
    plot(axM, M, errM(:,j), [ls 'o'], 'Color', color_of(pairRatios(j)), 'MarkerSize', 4, ...
        'MarkerFaceColor', P.bg, 'HandleVisibility', 'off');
end
for i = 1:numel(gaps)
    plot(axM, NaN, NaN, styles{i}, 'Color', P.fg, ...
        'DisplayName', sprintf('{\\it d} = %g', gaps(i)));
end

xlabel(axN, 'boundary resolution, {\it N_t}');
ylabel(axN, 'error, |{\it u}_{\rm cm}^{{\it N_t}} - {\it u}_{\rm cm}^*|');
axN.YScale = 'log';
xlim(axN, [0 110]); ylim(axN, [1e-16 1]); yticks(axN, 10.^(-16:4:0));

xlabel(axM, 'expansion order, {\it M}');
ylabel(axM, 'error, |{\it U}^{\it M} - {\it U}^{\rm ref}|_\infty');
axM.YScale = 'log';
xlim(axM, [0 max(M) + 5]); ylim(axM, [1e-14 1e-2]); yticks(axM, 10.^(-14:4:-2));

legend(axN, 'Location', 'northeast', 'Color', 'none');
legend(axM, 'Location', 'southwest', 'Color', 'none');

% ── Layout ────────────────────────────────────────────────────────────────────
ml = 1.25; mr = 0.25; mb = 0.9; mt = 0.3; mid = 1.5;

fig.Units = 'inches'; figW = fig.Position(3); figH = fig.Position(4);
axN.Units = 'inches'; axM.Units = 'inches';

panelW = (figW - ml - mr - mid) / 2;
panelH = figH - mb - mt;

axN.Position = [ml, mb, panelW, panelH];
axM.Position = [ml + panelW + mid, mb, panelW, panelH];

panel_label('\it (a)', axN);
panel_label('\it (b)', axM);

% ── Helpers ───────────────────────────────────────────────────────────────────
function [n, err, hdr] = load_errors(fname)
    lines = readlines(fname);
    lines = lines(startsWith(lines, '#'));
    hdr.ratios = header_values(lines, 'aspect ratios a/b:');
    hdr.gaps   = header_values(lines, 'edge gaps:');
    A   = readmatrix(fname, 'FileType', 'text', 'NumHeaderLines', numel(lines));
    n   = A(:,1);
    err = A(:,2:end);
end

function v = header_values(lines, key)
    tok = regexp(lines, ['^#\s*' key '\s*(.*)$'], 'tokens', 'once');
    tok = tok(~cellfun(@isempty, tok));
    v = [];
    if ~isempty(tok), v = str2num(tok{1}{1}); end %#ok<ST2NM>
end

function s = ratio_label(r)
    if r >= 1
        s = sprintf('%g', r);
    else
        s = sprintf('1/%g', 1/r);
    end
end

function a = panel_label(str, ax)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', ax.FontSize, 'Units', 'inches');
    a.Position(1) = ax.Position(1) - 1.2;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + 0.1;
end
