% fig_tandem.m — far-field tandem chains across aspect ratio, started from an orientation
% wave along the chain: kymographs of the orientations (top); growth rates and frequencies
% of the transverse modes of the infinite chain against mode number, from
% farfield_tandem_stability.jl, and the Fourier amplitude of the imposed mode over time
% (bottom).

close all

P = palette();

ab_vals = [0.25, 1.0, 4.0];
ab_strs = {'0.25', '1.0', '4.0'};
m       = 10;   % imposed mode
cols    = P.cmap(7);  cols = cols([3 5 7], :);
outDir  = fullfile(fileparts(mfilename('fullpath')), '..', 'scripts', 'farfield', 'output');
Nc = numel(ab_vals);

fig = journal_figure([6.75 4.2], 2, 'tex');

ax = gobjects(2, Nc);  cb = gobjects(1, Nc);  thmax = zeros(1, Nc);
t_max = 2000;
for k = 1:numel(ax), ax(k) = axes(fig); hold(ax(k), 'on'); box(ax(k), 'on'); end

for c = 1:Nc
    % ── Simulation ────────────────────────────────────────────────────────────
    dir_c  = fullfile(outDir, ['tandem_' ab_strs{c} '_N100_d10.0_m10_ε0.01']);
    f      = dir(fullfile(dir_c, 'frame_*.dat'));
    frames = cellfun(@(f) load(fullfile(dir_c, f)), sort({f.name}), 'UniformOutput', false);  % columns: t  x  y  θ  vx  vy  ω
    t      = cellfun(@(d) d(1,1), frames);
    th     = cell2mat(cellfun(@(d) d(:,4), frames, 'UniformOutput', false));  % particle × time
    N      = size(th, 1);
    A      = 2 * abs(fft(th, [], 1)) / N;   % Fourier amplitudes over the particle index

    imagesc(ax(1,c), t, 1:N, th);
    axis(ax(1,c), 'tight');
    colormap(ax(1,c), cmocean('phase'));
    thmax(c) = max(abs(th(:)));
    cb(c) = colorbar(ax(1,c), 'Location', 'northoutside');
    cb(c).Label.String = '{\it\theta}';
    set(ax(1,c), 'Layer', 'top');
    xlabel(ax(1,c), 'time, {\it t}');
    text(ax(1,c), 0.02, 0.97, sprintf('\\it (%s) {\\it a/b} \\rm = %.2f', char('a' + c - 1), ab_vals(c)), ...
         'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', ax(1,c).FontSize, ...
         'BackgroundColor', 'w', 'Margin', 1);

    name = sprintf('{\\it a/b} = %.2f', ab_vals(c));
    plot(ax(2,3), t, A(m+1,:), '-', 'Color', cols(c,:), 'DisplayName', name);
    if t(end) < 4000   % stopped on contact
        plot(ax(2,3), t(end), A(m+1,end), 'o', 'Color', cols(c,:), 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    end

    % ── Linear stability ──────────────────────────────────────────────────────
    % a real perturbation has the conjugate pairs Re λ± ± i |Im λ±|
    D = load(fullfile(outDir, ['tandem_stability_' ab_strs{c} '.dat']));  % columns: m  κ  Re λ₊  |Im λ₊|  Re λ₋  |Im λ₋|
    plot(ax(2,1), D(:,1), D(:,3), '-',  'Color', cols(c,:));
    plot(ax(2,1), D(:,1), D(:,5), '--', 'Color', cols(c,:));
    plot(ax(2,2), D(:,1), [D(:,4), -D(:,4)], '-',  'Color', cols(c,:));
    plot(ax(2,2), D(:,1), [D(:,6), -D(:,6)], '--', 'Color', cols(c,:));
end

% oblate and disk share the scale of the imposed wave; prolate grows to O(1)
for c = 1:2, clim(ax(1,c), max(thmax(1:2)) * [-1 1]); end
clim(ax(1,3), pi/2 * [-1 1]);
set(cb(3), 'Ticks', pi/2 * [-1 0 1], 'TickLabels', {'-\pi/2', '0', '\pi/2'});
xlim(ax(1,:), [0 t_max]);  xlim(ax(2,3), [0 t_max]);

ylabel(ax(1,1), 'particle, {\it i}');
set(ax(1,2:end), 'YTickLabel', {});
xlabel(ax(2,1), 'mode, {\it m}');  ylabel(ax(2,1), 'Re {\it\lambda}');
xlabel(ax(2,2), 'mode, {\it m}');  ylabel(ax(2,2), 'Im {\it\lambda}');
xlabel(ax(2,3), 'time, {\it t}');  ylabel(ax(2,3), sprintf('|{\\it\\theta}_{%d}|', m));
set(ax(2,3), 'YScale', 'log');
xlim(ax(2,1:2), [0 N/2]);
ylim(ax(2,2), 1.15 * max(abs(ylim(ax(2,2)))) * [-1 1]);
ylim(ax(2,3), [1e-5 1]);  yticks(ax(2,3), [1e-4 1e-2 1]);
legend(ax(2,3), 'Location', 'southeast', 'Color', 'none');
for c = 1:Nc
    text(ax(2,c), 0.02, 0.97, sprintf('\\it (%s)', char('a' + Nc + c - 1)), 'Units', 'normalized', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', ax(2,c).FontSize);
end

% ── Layout ────────────────────────────────────────────────────────────────────

figW = fig.PaperSize(1);  figH = fig.PaperSize(2);
lm = 1.3; rm = 0.4; bm = 0.9; tm = 1.1; cg = 1.3; rg = 1.0;
panelW = (figW - lm - rm - (Nc-1)*cg) / Nc;
panelH = (figH - bm - tm - rg) / 2;
set(ax, 'Units', 'inches', 'ActivePositionProperty', 'position');
set(cb, 'Units', 'inches');
for c = 1:Nc
    x0 = lm + (c-1)*(panelW + cg);
    ax(1,c).Position = [x0, bm + panelH + rg,              panelW, panelH];
    ax(2,c).Position = [x0, bm,                            panelW, panelH];
    cb(c).Position   = [x0, bm + 2*panelH + rg + 0.1,      panelW, 0.15];
end
