% meanSquareDisplacement.m — mean-squared displacement from saved suspension
% frames, checking ballistic (slope 2) scaling, alongside configuration
% snapshots fading from transparent (early) to opaque (late).

close all

frameDir = fullfile(fileparts(mfilename('fullpath')), '..', 'examples', 'suspension', 'frames_1.0');
files    = dir(fullfile(frameDir, 'frame_*.dat'));
names    = sort({files.name});
Nframes  = numel(names);

t = zeros(Nframes, 1);
X = [];  Y = [];
for k = 1:Nframes
    A    = readmatrix(fullfile(frameDir, names{k}));  % columns: t x y θ vx vy ω
    t(k) = A(1, 1);
    X    = [X, A(:,2)]; %#ok<AGROW>
    Y    = [Y, A(:,3)]; %#ok<AGROW>
end
X = X.';  Y = Y.';   % Nframes × Nparticles

x0  = X(1, :);
y0  = Y(1, :);
msd = mean((X - x0).^2 + (Y - y0).^2, 2);

cx = mean(X, 2);
cy = mean(Y, 2);
R  = max(hypot(X - cx, Y - cy), [], 2);

idx     = t > 0;
t_fit   = t(idx);
msd_fit = msd(idx);
R_fit   = R(idx);

radius     = 1.0;
Nsnaps     = 5;
snapIdx    = round(linspace(1, Nframes, Nsnaps));
snapAlphas = linspace(0.2, 1.0, Nsnaps);

cmap      = cmocean('tempo', 256);
snapCidx  = max(1, round(255 * t(snapIdx) / t(end))) + 1;
snapColor = cmap(snapCidx, :);

% --- Figure -------------------------------------------------------------
fig = journal_figure([6.75, 2.25], 2, 'tex');
P   = palette();
col = P.col;

axCfg = axes(fig); hold(axCfg, 'on'); box(axCfg, 'on'); axis(axCfg, 'equal');
axR   = axes(fig); hold(axR,   'on'); box(axR,   'on');
axMSD = axes(fig); hold(axMSD, 'on'); box(axMSD, 'on');

for n = 1:Nsnaps
    k = snapIdx(n);
    draw_config(axCfg, X(k,:), Y(k,:), radius, snapColor(n,:), snapAlphas(n));
end
xlabel(axCfg, '{\it x}'); ylabel(axCfg, '{\it y}');
axis(axCfg, 50 * [-1 1 -1 1]);

colormap(axCfg, cmap);
clim(axCfg, [0, t(end)]);
cb = colorbar(axCfg, 'horizontal');
cb.Label.String = 'time, {\it t}';

plot(axR, t_fit, R_fit/R_fit(1), '-', 'Color', col, 'MarkerFaceColor', col);
xlabel(axR, 'time, {\it t}'); ylabel(axR, 'radius, {\it R}/{\it R}_0');


plot(axMSD, t_fit, msd_fit, 'o-', 'Color', col, 'MarkerFaceColor', col, 'DisplayName', 'simulation');
plot(axMSD, t_fit, 0.1 * msd_fit(end) * (t_fit/t_fit(end)).^2, 'k--', 'DisplayName', '\propto{\it t}^2');
plot(axMSD, t_fit, 0.098175^2 * t_fit.^2, 'k--', 'DisplayName', '\propto{\it t}^2');

set(axMSD, 'XScale', 'log', 'YScale', 'log');
xlabel(axMSD, 'time, {\it t}'); ylabel(axMSD, 'MSD');
legend(axMSD, 'Location', 'southeast');


% --- Layout ---------------------------------------------------------------
lm = 0.6; rm = 0.4; bm = 1.0; tm = 0.20; gapH = 1;

fig.Units = 'inches'; figW = fig.Position(3); figH = fig.Position(4);
set([axCfg, axR, axMSD], 'Units', 'inches');
panelW = (figW - lm - rm - 2 * gapH) / 3;
panelH = figH - bm - tm;

axCfg.Position = [lm,                         bm, panelW, panelH];
axR.Position   = [lm + panelW + gapH,         bm, panelW, panelH];
axMSD.Position = [lm + 2 * (panelW + gapH),   bm, panelW, panelH];

cb.Units    = 'inches';
cb.Position = [lm + 0.15, bm + 0.30, 0.45 * panelW, 0.10];

panel_label('\it (a)', axCfg);
panel_label('\it (b)', axR);
panel_label('\it (c)', axMSD);

% --- helpers ------------------------------------------------------------
function a = panel_label(str, ax)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', 18, 'Units', 'inches');
    a.Position(1) = ax.Position(1) - 0.05;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4);
end

function slope_line(ax, t, y, slope)
    offset = 0.5;   % shift below the data so the two don't overlap
    y_ref  = offset * y(end) * (t/t(end)).^slope;
    plot(ax, t, y_ref, 'k--', 'HandleVisibility', 'off');
end

function draw_config(ax, x, y, r, faceColor, faceAlpha)
    phi = linspace(0, 2*pi, 50);
    for i = 1:numel(x)
        patch(ax, x(i) + r*cos(phi), y(i) + r*sin(phi), faceColor, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', 'none');
    end
end
