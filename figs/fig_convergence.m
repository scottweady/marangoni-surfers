% fig_convergence.m — convergence of the single-particle solver against the
% semi-analytical center-of-mass velocity, assessed independently with respect
% to the expansion order M and the boundary resolution N_t, for several particle
% aspect ratios a/b at fixed area (ab = 1).
%
% Panel (a): error vs expansion order M, boundary held at N_t = N_fine.
% Panel (b): error vs boundary resolution N_t, expansion held at M = M_fine.
%
% Requires: scripts/convergence/output/convergence_expansion.dat
%           scripts/convergence/output/convergence_boundary.dat

close all;

baseDir = fullfile(dataRoot, 'convergence');

[M, errM, ratios] = load_errors(fullfile(baseDir, 'convergence_expansion.dat'));
[N, errN, ~     ] = load_errors(fullfile(baseDir, 'convergence_boundary.dat'));

% ── Colors ────────────────────────────────────────────────────────────────────
nr   = numel(ratios);
P = palette();
cmap = P.cmap(nr + 2);
cmap = cmap(2:nr+1, :);

% ── Figure ────────────────────────────────────────────────────────────────────
fig = journal_figure([6.75, 2.25], 2, 'tex');

axM = axes(fig); hold(axM, 'on'); box(axM, 'on');
axN = axes(fig); hold(axN, 'on'); box(axN, 'on');

for k = 1:nr
    plot(axM, M, errM(:,k), '-o', 'Color', cmap(k,:), 'MarkerSize', 4, ...
        'MarkerFaceColor', 'w', 'DisplayName', ['{\it a/b} = ' ratio_label(ratios(k))]);
    plot(axN, N, errN(:,k), '-o', 'Color', cmap(k,:), 'MarkerSize', 4, ...
        'MarkerFaceColor', 'w', 'DisplayName', ['{\it a/b} = ' ratio_label(ratios(k))]);
end

xlabel(axM, 'expansion order, {\it M}');
ylabel(axM, 'error, |{\it u}_{\rm cm}^{\it M} - {\it u}_{\rm cm}^*|');
style_log_axis(axM, M, errM);

xlim(axM, [0 28])
ylim(axM, [1e-17 1e-3])
yticks(axM, 10.^(-16:4:0))

xlabel(axN, 'boundary resolution, {\it N_t}');
ylabel(axN, 'error, |{\it u}_{\rm cm}^{{\it N_t}} - {\it u}_{\rm cm}^*|');

style_log_axis(axN, N, errN);
xlim(axN, [0 110])
ylim(axN, [1e-17 1])
yticks(axN, 10.^(-16:4:0))

legend(axM, 'Location', 'northeast', 'Color', 'none');
% legend(axN, 'Location', 'southwest', 'Color', 'none');

% ── Layout ────────────────────────────────────────────────────────────────────
ml = 1.25; mr = 0.25; mb = 0.9; mt = 0.3; mid = 1.5;

fig.Units = 'inches'; figW = fig.Position(3); figH = fig.Position(4);
axM.Units = 'inches'; axN.Units = 'inches';

panelW = (figW - ml - mr - mid) / 2;
panelH = figH - mb - mt;

axM.Position = [ml, mb, panelW, panelH];
axN.Position = [ml + panelW + mid, mb, panelW, panelH];

panel_label('\it (a)', axM);
panel_label('\it (b)', axN);

% ── Helpers ───────────────────────────────────────────────────────────────────
function [n, err, ratios] = load_errors(fname)
    fid = fopen(fname);
    if fid == -1
        error('Cannot open %s — run scripts/convergence/convergence.jl first.', fname);
    end
    ratios = [];
    nhdr   = 0;
    while true
        pos  = ftell(fid);
        line = fgetl(fid);
        if ~ischar(line) || isempty(line) || line(1) ~= '#'
            fseek(fid, pos, 'bof');
            break
        end
        nhdr = nhdr + 1;
        tok  = regexp(line, '^#\s*aspect ratios a/b:\s*(.*)$', 'tokens', 'once');
        if ~isempty(tok), ratios = str2num(tok{1}); end %#ok<ST2NM>
    end
    fclose(fid);

    A   = readmatrix(fname, 'FileType', 'text', 'NumHeaderLines', nhdr);
    n   = A(:,1);
    err = A(:,2:end);
    if isempty(ratios), ratios = 1:size(err,2); end
end

function style_log_axis(ax, n, err)
    ax.YScale = 'log';
    xlim(ax, [0, max(n) + min(diff(n))]);

    lo = floor(log10(min(err(err > 0))));
    ylim(ax, [10^lo, 1]);

    step = 2 * ceil(-lo / 8);   % keep ~4-5 labelled decades
    ax.YTick = 10 .^ (lo : step : 0);
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
