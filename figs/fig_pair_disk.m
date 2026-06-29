% Three-panel figure:
%   Left  (full height): dist(Γ1,Γ2) vs time for circular pair (a=1)
%   Right top:           ucm vs dist(Γ1,Γ2)  (finite-differenced from trajectory)
%   Right bottom:        omega vs dist(Γ1,Γ2) (finite-differenced from trajectory)

close all;

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');

% ---- Load trajectory (circular, a=1) ----------------------------------------
trajFile = fullfile(dataDir, 'disk', 'a1.0000_b1.0000_d16.0.dat');
T   = readmatrix(trajFile);
t   = T(:,1);
x1  = T(:,2);  th1 = T(:,4);
gap = T(:,8);

% Finite-difference velocities from trajectory
ucm = gradient(x1,  t);               % dx1/dt  (x-velocity of particle 1)
omg = gradient(unwrap(th1), t);        % dθ1/dt  (angular velocity)

% ---- Figure -----------------------------------------------------------------
fig = journal_figure([6.5, 2.5], 2, 'tex');

axDist = axes(fig); hold(axDist, 'on'); box(axDist, 'on');
axU    = axes(fig); hold(axU,    'on'); box(axU,    'on');
axW    = axes(fig); hold(axW,    'on'); box(axW,    'on');

% ---- Plot -------------------------------------------------------------------
P   = palette();
col = P.col;

tf = find(gap < 1e-2, 1);   % index of first point where gap < 0.01
plot(axDist, t(1:tf), gap(1:tf), '-', 'Color', col, 'DisplayName', 'simulation');

% Fit gap = gap(1)*(1 - t/ts)^0.5 for ts using first quarter of data
nq    = round(tf/4);
model = @(ts, tt) 14 * (1 - tt/ts).^(4/9);
ts0   = t(tf) * 1.1;
ts_fit = lsqcurvefit(model, ts0, t(1:nq), gap(1:nq), t(tf), [], ...
            optimoptions('lsqcurvefit','Display','off'));
t_ext = linspace(t(1), ts_fit, 500);
plot(axDist, t_ext, model(ts_fit, t_ext), '--', 'Color', col, ...
    'DisplayName', '{\it d}_0(1-{\it t}/{\it t}_s)^{4/9}');
legend(axDist, 'Location', 'northeast');

xlabel(axDist, 'time, {\it t}');
ylabel(axDist, 'dist({\it\Gamma}_1, {\it\Gamma}_2)');
xlim(axDist, [0 9e3])

plot(axU, gap, ucm, '-', 'Color', col, 'DisplayName', 'simulation');
plot(axU, gap, ucm(1) * (gap(1)./gap).^(5/4), '--', 'Color', col, ...
    'DisplayName', '\propto{\it r}_{12}^{-5/4}');
legend(axU, 'Location', 'northeast');
ylabel(axU, '{\it u}_{cm}');
axU.XTickLabel = {};
xlim(axU, [0 1e1]);
ylim(axU, [0 0.01]);

plot(axW, gap, omg, '-', 'Color', col, 'DisplayName', 'simulation');
plot(axW, gap, 0*gap, ':', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off');
ylabel(axW, '{\it \omega}_{cm}');
xlabel(axW, 'dist({\it\Gamma}_1, {\it\Gamma}_2)');
% axU.XScale = 'log';
% axW.XScale = 'log';
linkaxes([axU, axW], 'x');
xlim(axW, [0 1e1]);
ylim(axW, [-2e-4 6e-4]);
yticks(axW, -2e-4:2e-4:6e-4);

% ---- Layout -----------------------------------------------------------------
for ax = [axDist, axU, axW]
    ax.ActivePositionProperty = 'position';
    ax.Units = 'inches';
end
drawnow;

tiD = axDist.TightInset;
tiU = axU.TightInset;
tiW = axW.TightInset;

lm   = max([0.8,  tiD(1), tiU(1), tiW(1)]) - 0.3;
rm   = max([0.15, tiD(3), tiU(3), tiW(3)]) + 0.06;
bm   = max([0.5,  tiD(2), tiW(2)]) - 0.05;
tm   = max([0.15, tiD(4), tiU(4)]) + 0.06;
midV = 0.5;   % vertical gap between right panels
midH = 1.2;   % horizontal gap between left and right

figW = fig.PaperSize(1);
figH = fig.PaperSize(2);

panelH_full = figH - bm - tm;
panelW_L    = 0.45 * (figW - lm - rm - midH);
panelW_R    = figW - lm - rm - midH - panelW_L;
panelH_R    = (panelH_full - midV) / 2;

xL = lm;
xR = lm + panelW_L + midH;

axDist.Position = [xL, bm,                         panelW_L, panelH_full];
axW.Position    = [xR, bm,                         panelW_R, panelH_R];
axU.Position    = [xR, bm + panelH_R + midV,       panelW_R, panelH_R];

drawnow;
axDist.Position = [xL, bm,                         panelW_L, panelH_full];
axW.Position    = [xR, bm,                         panelW_R, panelH_R];
axU.Position    = [xR, bm + panelH_R + midV,       panelW_R, panelH_R];

fs = axDist.FontSize;
axU.FontSize = fs;
axW.FontSize = fs;

% ---- Inset schematic: two disks ---------------------------------------------
r_disk  = 1.0;
gap_vis = 0.5;
xc      = [-r_disk - gap_vis/2,  r_disk + gap_vis/2];
phi_c   = linspace(0, 2*pi, 200);
alp     = pi/2;
phi_act = linspace(alp, 2*pi - alp, 150);

theta_disk = pi/2;
ct = cos(theta_disk); st = sin(theta_disk);

axIn = axes(fig, 'Units', 'inches');
hold(axIn, 'on'); axis(axIn, 'equal'); axis(axIn, 'off');
for k = 1:2
    fill(axIn, r_disk*cos(phi_c) + xc(k), r_disk*sin(phi_c), ...
        [0.92 0.92 0.92], 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 2);
    % Active arc rotated by theta=pi/2
    xa = r_disk*cos(phi_act); ya = r_disk*sin(phi_act);
    plot(axIn, ct*xa - st*ya + xc(k), st*xa + ct*ya, ...
        '-', 'Color', col, 'LineWidth', 4);
end
axIn.XLim = [xc(1) - r_disk - 0.2,  xc(2) + r_disk + 0.2];
axIn.YLim = [-r_disk - 0.2,  r_disk + 0.2];

% Position inset in upper-right corner of panel (a) after layout
insetW = 3;
insetH = insetW * diff(axIn.YLim) / diff(axIn.XLim);
insetX = axDist.Position(1) + 0.05;
insetY = axDist.Position(2) + 0.05;
axIn.Position = [insetX, insetY, insetW, insetH];

% Panel labels.
panel_label('\it (a)', axDist);
panel_label('\it (b)', axU);
panel_label('\it (c)', axW);

% ---- helpers ----------------------------------------------------------------
function a = panel_label(str, ax)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', 18, 'Units', 'inches');
    a.Position(1) = ax.Position(1) - 0.05;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4);
end
