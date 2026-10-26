% fig_patches_farfield.m — expanding circular patches across aspect ratio from direct
% simulation (top) and the far-field model (bottom), both started from the same
% configuration: initial and expanded configurations with particle trajectories.

close all

P = palette();
snap_cols = [P.snap(1,:); P.snap(3,:)];

ab_vals = [0.25, 1.0, 4.0];
ab_strs = {'0.25', '1.0', '4.0'};
models  = {'frames', 'farfield'};
labels  = {'direct', 'far-field'};
Nc = numel(ab_vals);  Nr = numel(models);

fig = journal_figure([6.75 4.5], 2, 'tex');

% ── Load frames ───────────────────────────────────────────────────────────────
% Each row is shown at its earliest final time across aspect ratios, the collision time of
% the prolate case.
allf = cell(Nr, Nc);  tv = cell(Nr, Nc);
for r = 1:Nr
    for c = 1:Nc
        dir_rc    = fullfile(dataRoot, 'patch', [models{r} '_' ab_strs{c}]);
        f         = dir(fullfile(dir_rc, 'frame_*.dat'));
        allf{r,c} = cellfun(@(f) load(fullfile(dir_rc, f)), sort({f.name}), 'UniformOutput', false);  % columns: t  x  y  θ  vx  vy  ω
        tv{r,c}   = cellfun(@(d) d(1,1), allf{r,c});
    end
end

D0 = cell(Nr, Nc);  DN = cell(Nr, Nc);  X = cell(Nr, Nc);  Y = cell(Nr, Nc);  t_show = zeros(Nr, Nc);
for r = 1:Nr
    t_snap = min(cellfun(@(t) t(end), tv(r,:)));
    for c = 1:Nc
        [~, kN] = min(abs(tv{r,c} - t_snap));
        frames  = allf{r,c}(1:kN);
        t_show(r,c) = tv{r,c}(kN);
        D0{r,c} = frames{1};
        DN{r,c} = frames{end};
        X{r,c}  = cell2mat(cellfun(@(d) d(:,2), frames, 'UniformOutput', false));  % particle × time
        Y{r,c}  = cell2mat(cellfun(@(d) d(:,3), frames, 'UniformOutput', false));
    end
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
ax = gobjects(Nr, Nc);

for r = 1:Nr
    for c = 1:Nc
        a = sqrt(ab_vals(c));  b = 1/a;
        ex = a*cos(s);  ey = b*sin(s);
        d0 = D0{r,c};  dN = DN{r,c};

        ax(r,c) = axes(fig); hold(ax(r,c), 'on'); axis(ax(r,c), 'equal');
        box(ax(r,c), 'on'); set(ax(r,c), 'XTick', [], 'YTick', []);

        % shaded extents: initial disk, expanded annulus
        xbar0 = mean(d0(:,2));  ybar0 = mean(d0(:,3));
        xbarN = mean(dN(:,2));  ybarN = mean(dN(:,3));
        R0_max = max(hypot(d0(:,2)-xbar0, d0(:,3)-ybar0));
        RN_min = min(hypot(dN(:,2)-xbarN, dN(:,3)-ybarN));
        RN_max = max(hypot(dN(:,2)-xbarN, dN(:,3)-ybarN));

        fill(ax(r,c), xbar0 + R0_max*cos(th_ring), ybar0 + R0_max*sin(th_ring), ...
             snap_cols(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        for sgn = [1, -1]
            th_half = linspace(0, sgn*pi, 128);
            fill(ax(r,c), [xbarN + RN_max*cos(th_half), xbarN + RN_min*cos(fliplr(th_half))], ...
                          [ybarN + RN_max*sin(th_half), ybarN + RN_min*sin(fliplr(th_half))], ...
                 snap_cols(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        end

        plot(ax(r,c), X{r,c}', Y{r,c}', '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.75);

        frames = {d0, dN};
        for f = 1:2
            fdata = frames{f};
            for i = 1:size(fdata, 1)
                R  = [cos(fdata(i,4)) -sin(fdata(i,4)); sin(fdata(i,4)) cos(fdata(i,4))];
                xy = R * [ex; ey];
                fill(ax(r,c), fdata(i,2)+xy(1,:), fdata(i,3)+xy(2,:), snap_cols(f,:), 'EdgeColor', 'none');
            end
        end

        xlim(ax(r,c), xlims + [-xpad, xpad]); ylim(ax(r,c), ylims + [-ypad, ypad]);

        text(ax(r,c), 0, 0.92, sprintf('\\it (%s) {\\it a/b} \\rm = %.2f, %s', char('a' + (r-1)*Nc + c - 1), ab_vals(c), labels{r}), ...
             'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', ax(r,c).FontSize);
    end
end

% scale bar in (a)
sb_len = 25;
line(ax(1,1), xlims(1) + [0 sb_len], ylims(1)*[1 1], 'Color', 'k', 'LineWidth', 2);
text(ax(1,1), xlims(1) + sb_len/2, ylims(1) + 0.02*diff(ylims), '25{\it{l}}_*', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', ax(1,1).FontSize);

% ── Layout ────────────────────────────────────────────────────────────────────
figW = fig.PaperSize(1);  figH = fig.PaperSize(2);
lm = 0.1; rm = 0.1; bm = 0.1; tm = 0.2; cg = 0.2; rg = 0.3;
panelW = (figW - lm - rm - (Nc-1)*cg) / Nc;
panelH = (figH - bm - tm - (Nr-1)*rg) / Nr;
set(ax, 'Units', 'inches', 'ActivePositionProperty', 'position');
for r = 1:Nr
    for c = 1:Nc
        ax(r,c).Position = [lm + (c-1)*(panelW + cg), bm + (Nr-r)*(panelH + rg), panelW, panelH];
    end
end
