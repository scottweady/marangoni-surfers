% fig_tandem.m — tandem chains across aspect ratio, started from an orientation
% wave along the chain: snapshots of the chains on the surfactant concentration, from
% tandem_snapshot.jl, in the frame of the chain centroid (top); growth rates and
% frequencies of the transverse modes of the infinite chain against mode number, from
% farfield_tandem_stability.jl, and the Fourier amplitude of the imposed mode over time
% against its linear prediction (bottom).

close all

P = palette();

ab_vals = [0.25, 1.0, 4.0];
ab_strs = {'0.25', '1.0', '4.0'};
m       = 5;   % imposed mode
cols    = P.cmap(7);  cols = cols([3 5 7], :);
t_snap  = [0 1000];   % snapshot times
outDir  = fullfile(dataRoot, 'tandem');
Nc = numel(ab_vals);  Ns = numel(t_snap);

% ── Snapshot grids ────────────────────────────────────────────────────────────
S = cell(Nc, Ns);
for c = 1:Nc
    for k = 1:Ns
        fid = fopen(fullfile(outDir, sprintf('tandem_snapshot_%s_t%d.dat', ab_strs{c}, t_snap(k))));
        hdr = fgetl(fid);
        Nx  = str2double(regexp(hdr, 'Nx=(\d+)', 'tokens', 'once'));
        Ny  = str2double(regexp(hdr, 'Ny=(\d+)', 'tokens', 'once'));
        S{c,k}.x   = str2num(fgetl(fid));  %#ok<ST2NM>
        S{c,k}.y   = str2num(fgetl(fid));  %#ok<ST2NM>
        S{c,k}.phi = fscanf(fid, '%f', [Nx, Ny])';
        fclose(fid);
    end
end

% ── Figure size ───────────────────────────────────────────────────────────────
% the snapshots are to scale, so their heights set the height of the figure
sc = 2;  figW = 6.75 * sc;
lm = 1.0; rm = 0.6; bm = 0.9; cg = 1; panelH = 2.7;
sg = 0.3; sr = 0.15; cbg = 0.2; cbw = 0.15; stm = 0.6; sb = 1;   % snapshots: column and row gaps, colorbar gap and width, top margin, gap below
snapW = (figW - lm - rm - (Ns-1)*sg - cbg - cbw) / Ns;   % spans the bottom row, with the colorbar
s     = snapW / (S{1,1}.x(end) - S{1,1}.x(1));          % inches per unit length
snapH = cellfun(@(Sc) (Sc.y(end) - Sc.y(1)) * s, S(:,1))';
figH  = bm + panelH + sb + sum(snapH) + (Nc-1)*sr + stm;

fig = journal_figure([figW figH] / sc, sc, 'tex');

ax = gobjects(1, Nc);  as = gobjects(Nc, Ns);
t_max = 500;
for k = 1:numel(ax), ax(k) = axes(fig); hold(ax(k), 'on'); box(ax(k), 'on'); end
for k = 1:numel(as), as(k) = axes(fig); hold(as(k), 'on'); box(as(k), 'on'); end
phi = linspace(0, 2*pi, 41);

% concentration colormap over [0, 1], varying only over c_rng and solid outside it
c_rng = [0.4 0.6];
cmap  = cmocean('dense', 256);
cmap  = cmap(round(1 + 255 * min(max((linspace(0, 1, 1001) - c_rng(1)) / diff(c_rng), 0), 1)), :);

for c = 1:Nc
    % ── Simulation ────────────────────────────────────────────────────────────
    dir_c  = fullfile(outDir, ['tandem_' ab_strs{c} '_N50_d10.0_m5_ε0.01']);
    f      = dir(fullfile(dir_c, 'frame_*.dat'));
    frames = cellfun(@(f) load(fullfile(dir_c, f)), sort({f.name}), 'UniformOutput', false);  % columns: t  x  y  θ  vx  vy  ω
    t      = cellfun(@(d) d(1,1), frames);
    th     = cell2mat(cellfun(@(d) d(:,4), frames, 'UniformOutput', false));  % particle × time
    N      = size(th, 1);
    A      = 2 * abs(fft(th, [], 1)) / N;   % Fourier amplitudes over the particle index

    % ── Snapshots ─────────────────────────────────────────────────────────────
    % the label sits above the chain, which starts on the centerline
    a = sqrt(ab_vals(c));  b = 1/a;
    for k = 1:Ns
        [~, n] = min(abs(t - t_snap(k)));
        d = frames{n};
        imagesc(as(c,k), S{c,k}.x, S{c,k}.y, S{c,k}.phi, 'AlphaData', ~isnan(S{c,k}.phi));
        for i = 1:N
            R  = [cos(d(i,4)) -sin(d(i,4)); sin(d(i,4)) cos(d(i,4))];
            xy = R * [a*cos(phi); b*sin(phi)];
            fill(as(c,k), d(i,2) - mean(d(:,2)) + xy(1,:), d(i,3) - mean(d(:,3)) + xy(2,:), 'w', ...
                 'EdgeColor', 'k', 'LineWidth', 0.25);
        end
        set(as(c,k), 'YDir', 'normal', 'DataAspectRatio', [1 1 1], 'XLim', S{c,k}.x([1 end]), ...
            'YLim', S{c,k}.y([1 end]), 'YTick', [], 'Layer', 'top');
        colormap(as(c,k), cmap);
        clim(as(c,k), [0 1]);
    end
    text(as(c,1), 0.01, 0.95, sprintf('{\\it a/b} = %.2f', ab_vals(c)), 'Units', 'normalized', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', as(c,1).FontSize, 'Color', 'w');

    name = sprintf('{\\it a/b} = %.2f', ab_vals(c));
    plot(ax(3), t, A(m+1,:), '-', 'Color', cols(c,:), 'DisplayName', name);
    if t(end) < t_snap(end)   % stopped on contact
        plot(ax(3), t(end), A(m+1,end), 'o', 'Color', cols(c,:), 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    end
    Q = load(fullfile(outDir, ['tandem_prediction_' ab_strs{c} '.dat']));  % columns: t  |Θ|
    plot(ax(3), Q(:,1), Q(:,2), ':', 'Color', cols(c,:), 'HandleVisibility', 'off');

    % ── Linear stability ──────────────────────────────────────────────────────
    % a real perturbation has the conjugate pairs Re λ± ± i |Im λ±|
    D = load(fullfile(outDir, ['tandem_stability_' ab_strs{c} '.dat']));  % columns: m  κ  Re λ₊  |Im λ₊|  Re λ₋  |Im λ₋|
    plot(ax(1), D(:,1), D(:,3), '-',  'Color', cols(c,:));
    plot(ax(1), D(:,1), D(:,5), '--', 'Color', cols(c,:));
    plot(ax(2), D(:,1), [D(:,4), -D(:,4)], '-',  'Color', cols(c,:));
    plot(ax(2), D(:,1), [D(:,6), -D(:,6)], '--', 'Color', cols(c,:));
end

plot(ax(3), NaN, NaN, 'k:', 'DisplayName', 'linear theory');

for k = 1:Ns, title(as(1,k), sprintf('{\\it t} = %d', t_snap(k)), 'FontWeight', 'normal'); end
set(as(1:end-1,:), 'XTickLabel', {});
for k = 1:Ns, xlabel(as(end,k), ['{\itx} - {\itx}{\it', char(772), '}']); end
text(as(1,1), 0, 1.05, '\it (a)', 'Units', 'normalized', 'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'bottom', 'FontSize', as(1,1).FontSize);
cbs = colorbar(as(end,end), 'Ticks', [0 c_rng 1]);
cbs.Title.String = '{\it c}';

xlim(ax(3), [0 t_max]);

xlabel(ax(1), 'mode, {\it m}');  ylabel(ax(1), 'Re {\it\lambda}');
xlabel(ax(2), 'mode, {\it m}');  ylabel(ax(2), 'Im {\it\lambda}');
xlabel(ax(3), 'time, {\it t}');  ylabel(ax(3), sprintf('|{\\it\\theta}_{%d}|', m));
xticks(ax(1:2), 0:5:25)
set([ax(1).YAxis ax(2).YAxis], 'Exponent', -3);   % short tick labels, to keep the panels close

set(ax(3), 'YScale', 'log');
xlim(ax(1:2), [0 N/2]);
ylim(ax(2), 1.15 * max(abs(ylim(ax(2)))) * [-1 1]);
ylim(ax(3), [1e-5 1]);  yticks(ax(3), [1e-4 1e-2 1]);
legend(ax(3), 'Location', 'southeast', 'Color', 'none');
for c = 1:Nc
    text(ax(c), 0, 0.99, sprintf('\\it (%s)', char('a' + c)), 'Units', 'normalized', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', ax(c).FontSize);
end

% ── Layout ────────────────────────────────────────────────────────────────────

panelW = (figW - lm - rm - (Nc-1)*cg) / Nc;
set([ax(:); as(:)], 'Units', 'inches', 'ActivePositionProperty', 'position');
set(cbs, 'Units', 'inches');
for c = 1:Nc
    ax(c).Position = [lm + (c-1)*(panelW + cg), bm, panelW, panelH];
end
y0 = figH - stm;
for c = 1:Nc
    y0 = y0 - snapH(c);
    for k = 1:Ns
        as(c,k).Position = [lm + (k-1)*(snapW + sg), y0, snapW, snapH(c)];
    end
    y0 = y0 - sr;
end
cbs.Position = [figW - rm - cbw, y0 + sr, cbw, figH - stm - y0 - sr];
