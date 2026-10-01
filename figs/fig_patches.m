% fig_suspension_patches.m — expanding circular patches across aspect ratio: initial
% (t = 0) and expanded (t = t_snap) particle configurations, in the style of
% fig_expansion panel (a).

close all

baseDir = fileparts(mfilename('fullpath'));
P = palette();
snap_cols = [P.snap(1,:); P.snap(3,:)];

ab_vals = [0.25, 1.0, 4.0];
ab_dirs = {'frames_0.25', 'frames_1.0', 'frames_4.0'};
Nc      = numel(ab_vals);

fig = journal_figure([6.75 2.4], 2, 'tex');

% ── Load frames up to the earliest final time across runs ─────────────────────
dataDirs = fullfile(baseDir, '..', 'data', 'patch', ab_dirs);
files = cell(Nc, 1);  tv = cell(Nc, 1);
for c = 1:Nc
    f        = dir(fullfile(dataDirs{c}, 'frame_*.dat'));
    files{c} = sort({f.name});
    tv{c}    = cellfun(@(f) dlmread(fullfile(dataDirs{c}, f), ' ', [0 0 0 0]), files{c});
end
t_snap = min(cellfun(@(t) t(end), tv));

D0 = cell(Nc, 1);  DN = cell(Nc, 1);  X = cell(Nc, 1);  Y = cell(Nc, 1);
for c = 1:Nc
    [~, kN] = min(abs(tv{c} - t_snap));
    frames  = cellfun(@(f) load(fullfile(dataDirs{c}, f)), files{c}(1:kN), 'UniformOutput', false);  % columns: t  x  y  θ  vx  vy  ω
    D0{c}   = frames{1};
    DN{c}   = frames{end};
    X{c}    = cell2mat(cellfun(@(d) d(:,2), frames, 'UniformOutput', false));  % particle × time
    Y{c}    = cell2mat(cellfun(@(d) d(:,3), frames, 'UniformOutput', false));
end

% shared limits from the largest expanded patch
allN  = vertcat(DN{:});
pad   = 2 * max(sqrt(ab_vals));
xlims = [min(allN(:,2))-pad, max(allN(:,2))+pad];
ylims = [min(allN(:,3))-pad, max(allN(:,3))+pad];
xpad  = 0.05 * diff(xlims);  ypad = 0.05 * diff(ylims);

% ── Panels ────────────────────────────────────────────────────────────────────
s  = linspace(0, 2*pi, 64);
th_ring = linspace(0, 2*pi, 256);
ax = gobjects(Nc, 1);

for c = 1:Nc
    a = sqrt(ab_vals(c));  b = 1/a;
    ex = a*cos(s);  ey = b*sin(s);
    d0 = D0{c};  dN = DN{c};

    ax(c) = axes(fig); hold(ax(c), 'on'); axis(ax(c), 'equal');
    box(ax(c), 'on'); set(ax(c), 'XTick', [], 'YTick', []);

    % shaded extents: initial disk, expanded annulus
    xbar0 = mean(d0(:,2));  ybar0 = mean(d0(:,3));
    xbarN = mean(dN(:,2));  ybarN = mean(dN(:,3));
    R0_max = max(hypot(d0(:,2)-xbar0, d0(:,3)-ybar0));
    RN_min = min(hypot(dN(:,2)-xbarN, dN(:,3)-ybarN));
    RN_max = max(hypot(dN(:,2)-xbarN, dN(:,3)-ybarN));

    fill(ax(c), xbar0 + R0_max*cos(th_ring), ybar0 + R0_max*sin(th_ring), ...
         snap_cols(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    for sgn = [1, -1]
        th_half = linspace(0, sgn*pi, 128);
        fill(ax(c), [xbarN + RN_max*cos(th_half), xbarN + RN_min*cos(fliplr(th_half))], ...
                    [ybarN + RN_max*sin(th_half), ybarN + RN_min*sin(fliplr(th_half))], ...
             snap_cols(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    end

    % trajectories
    plot(ax(c), X{c}', Y{c}', '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.75);

    % particles
    frames = {d0, dN};
    for f = 1:2
        fdata = frames{f};
        for i = 1:size(fdata, 1)
            R  = [cos(fdata(i,4)) -sin(fdata(i,4)); sin(fdata(i,4)) cos(fdata(i,4))];
            xy = R * [ex; ey];
            fill(ax(c), fdata(i,2)+xy(1,:), fdata(i,3)+xy(2,:), snap_cols(f,:), 'EdgeColor', 'none');
        end
    end

    xlim(ax(c), xlims + [-xpad, xpad]); ylim(ax(c), ylims + [-ypad, ypad]);

    text(ax(c), 0.03, 0.97, sprintf('\\it (%s) {\\it a/b} \\rm = %.2f', char('a' + c - 1), ab_vals(c)), ...
         'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
         'FontSize', ax(c).FontSize);
end

% scale bar in (a)
sb_len = 25;
line(ax(1), xlims(1) + [0 sb_len], ylims(1)*[1 1], 'Color', 'k', 'LineWidth', 2);
text(ax(1), xlims(1) + sb_len/2, ylims(1) + 0.02*diff(ylims), '25{\it R}', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', ax(1).FontSize);

% ── Layout ────────────────────────────────────────────────────────────────────
figW = fig.PaperSize(1);  figH = fig.PaperSize(2);
lm = 0.1; rm = 0.1; bm = 0.1; tm = 0.1; cg = 0.2;
panelW = (figW - lm - rm - (Nc-1)*cg) / Nc;
panelH = figH - bm - tm;
set(ax, 'Units', 'inches', 'ActivePositionProperty', 'position');
for c = 1:Nc
    ax(c).Position = [lm + (c-1)*(panelW + cg), bm, panelW, panelH];
end
