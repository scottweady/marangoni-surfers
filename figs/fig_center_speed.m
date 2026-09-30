% Two-panel figure (parallel configuration only):
%   Panel (a): concentration field for theta=pi/2  (N=7, d=10)
%   Panel (b): center speed vs N for d = 5, 10, 15, 20 with multipole prediction
%
% Requires: center_speed_snapshot_thpi2.dat, center_speed_vs_Nd.dat

close all;

baseDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'collective');

% ---- Speed data -------------------------------------------------------------
A    = readmatrix(fullfile(baseDir, 'center_speed_vs_Nd.dat'), 'NumHeaderLines', 1);
N    = A(:,1);
uN   = A(:,2:end);
dvals = [5 10 15 20];

uS   = 0.098175;
beta = 1.570796;
u_mp = zeros(numel(N), numel(dvals));
for j = 1:numel(dvals)
    u_mp(:,j) = arrayfun(@(Nk) mp_speed(Nk, dvals(j), uS, beta), N);
end

d    = 10;   % spacing in the snapshot

% ---- Colors -----------------------------------------------------------------
if ~exist('P', 'var'), P = palette(); end

% ---- Particle geometry (for overlays) ---------------------------------------
Np    = 7;
r     = 1.0;
xpos  = d * ((1:Np) - (Np+1)/2);
phi_c = linspace(0, 2*pi, 16);

% ==========================================================================
% Figure
% ==========================================================================
fig = journal_figure([6.5, 2.5], 2, 'tex');

axConc  = axes(fig); hold(axConc, 'on');
axSpeed = axes(fig); hold(axSpeed, 'on'); box(axSpeed, 'on');

% ==========================================================================
% Panel (a): concentration field
% ==========================================================================
[xplt, yplt, phi, U, V] = load_snapshot(fullfile(baseDir, 'center_speed_snapshot_thpi2.dat'));

imagesc(axConc, xplt, yplt, phi);
axConc.YDir = 'normal';
colormap(axConc, P.cmap());
clim(axConc, [0 1]);
axis(axConc, 'equal');
axConc.XTick = []; axConc.YTick = [];
box(axConc, 'on'); axConc.XColor = 'k'; axConc.YColor = 'k';

% Circle overlays
for k = 1:Np
    xe = r*cos(phi_c) + xpos(k);
    ye = r*sin(phi_c);
    patch(axConc, xe, ye, 'w', 'FaceColor', 'w', 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 2);
end

% Streamlines: seeds equispaced on x-axis (y=0), forward + backward
Nsl = 51;
sx  = linspace(xplt(1), xplt(end), Nsl);
sy  = zeros(1, Nsl);
[X, Y] = meshgrid(xplt, yplt);
h1 = streamline(axConc, X, Y,  U,  V, sx, sy, [1 100]);
h2 = streamline(axConc, X, Y, -U, -V, sx, sy, [1 100]);
set([h1; h2], 'Color', 'w', 'LineWidth', 1);

cb = colorbar(axConc, 'Location', 'westoutside', 'Ticks', 0:0.5:1);

% ==========================================================================
% Panel (b): center speed vs N
% ==========================================================================
plot(axSpeed, [0 14], [1 1], 'k:', 'HandleVisibility', 'off')
for j = 1:numel(dvals)
    col = P.lines(j,:);
    plot(axSpeed, N, uN(:,j), 'o-', 'Color', col, 'DisplayName', sprintf('{\\it d} = %d', dvals(j)));
    plot(axSpeed, N, u_mp(:,j), '--', 'Color', col, 'HandleVisibility', 'off');
end
plot(axSpeed, NaN, NaN, '--', 'Color', [0.3 0.3 0.3], 'DisplayName', 'multipole');

xlabel(axSpeed, 'number of particles, {\it{N}}');
ylabel(axSpeed, 'center speed, {\it{u}}_{cm} / {\it{U}}');
axSpeed.XTick = N';
legend(axSpeed, 'Location', 'northwest', 'Color', 'none');
xlim(axSpeed, [0 14]);
ylim(axSpeed, [0.9 2.5]);
% ==========================================================================
% Layout (same as fig_center_speed; concentration panel fills both strips)
% ==========================================================================
for ax = [axConc, axSpeed]
    ax.ActivePositionProperty = 'position';
    ax.Units = 'inches';
end
drawnow;

tiS = axSpeed.TightInset;
lm     = 0.8;
rm     = max(0.15, tiS(3)) + 0.06;
bm     = max(0.42, tiS(2)) + 0.08;
tm     = max(0.15, tiS(4)) + 0.06;
midH   = 1;
midV   = 0.25;

figW = fig.PaperSize(1);
figH = fig.PaperSize(2);

panelH_full = figH - bm - tm;

panelW_speed = 0.35 * (figW - lm - rm - midH);
panelW_conc  = figW - lm - rm - midH - panelW_speed;
panelH_conc  = 2 * (panelW_conc / 4) + midV;

xConc  = lm;
xSpeed = lm + panelW_conc + midH;

% Full x extent (all particles); y fit to the panel aspect ratio
xlim(axConc, [xplt(1) xplt(end)]);
ylim(axConc, (xplt(end) - xplt(1)) * panelH_conc / panelW_conc / 2 * [-1 1]);

axConc.Position  = [xConc,  bm, panelW_conc,  panelH_conc];
axSpeed.Position = [xSpeed, bm, panelW_speed, panelH_full];

drawnow;
axConc.Position  = [xConc,  bm, panelW_conc,  panelH_conc];
axSpeed.Position = [xSpeed, bm, panelW_speed, panelH_full];

cb.Units = 'inches';
cbW = 0.1;
cbG = 0.08;
cb.Position = [xConc - cbG - cbW, bm, cbW, panelH_conc];

cb.Title.String   = '{\it \phi}';
cb.Title.FontName = 'Times';
cb.Title.FontSize = axConc.FontSize;

panel_label('\it (a)', axConc);
panel_label('\it (b)', axSpeed);

% ==========================================================================
% Helpers
% ==========================================================================
function [xplt, yplt, phi, U, V] = load_snapshot(fname)
    fid = fopen(fname);
    if fid == -1
        error('Cannot open %s — run center_speed_snapshot.jl first.', fname);
    end
    hdr  = fgetl(fid);
    Nx   = str2double(regexp(hdr, 'Nx=(\d+)', 'tokens', 'once'));
    Ny   = str2double(regexp(hdr, 'Ny=(\d+)', 'tokens', 'once'));
    xplt = str2num(fgetl(fid)); %#ok<ST2NM>
    yplt = str2num(fgetl(fid)); %#ok<ST2NM>
    phi  = fscanf(fid, '%f', [Nx, Ny])';
    U    = fscanf(fid, '%f', [Nx, Ny])';
    V    = fscanf(fid, '%f', [Nx, Ny])';
    fclose(fid);
end

function u = mp_speed(N, d, uS, beta)
    ic  = (N + 1) / 2;
    xj  = d * ((1:N) - (N+1)/2);
    sep = abs(xj(ic) - xj([1:ic-1, ic+1:end]));
    u   = (uS + sum(beta ./ (4*pi*sep))) / uS;
end

function a = panel_label(str, ax)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', 20, 'Units', 'inches');
    a.Position(1) = ax.Position(1) - 0.9;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + 0.1;
end
