% Three-panel figure:
%   Panel (a): concentration field for theta=pi/2  (N=13, d=20)
%   Panel (b): concentration field for theta=0
%   Panel (c): center speed vs N with multipole prediction
%
% Requires: center_speed_snapshot_thpi2.dat, center_speed_snapshot_th0.dat
%           center_speed_vs_N_pi2.dat, center_speed_vs_N_0.dat

close all;

baseDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'collective');

% ---- Speed data -------------------------------------------------------------
A2  = readmatrix(fullfile(baseDir, 'center_speed_vs_N_pi2.dat'), 'NumHeaderLines', 1);
A0  = readmatrix(fullfile(baseDir, 'center_speed_vs_N_0.dat'),   'NumHeaderLines', 1);
N   = A2(:,1);
uN2 = A2(:,2);
uN0 = A0(:,2);

uS   = 0.098175;
beta = 1.570796;
d    = 20;
u_mp = arrayfun(@(Nk) mp_speed(Nk, d, uS, beta), N);

% ---- Colors -----------------------------------------------------------------
if ~exist('P', 'var'), P = palette(); end
c0 = P.c1;
c2 = P.c2;

% ---- Particle geometry (for overlays) ---------------------------------------
Np      = 13;
r       = 1.0;
xpos    = d * ((1:Np) - (Np+1)/2);
ic      = (Np+1)/2;
phi_c   = linspace(0, 2*pi, 16);
alp     = pi/2;
phi_act = linspace(alp, 2*pi-alp, 16);

% ==========================================================================
% Figure  — taller to stack two concentration strips
% ==========================================================================
% Pre-compute figure height: two 4:1 concentration strips stacked
% rough margins: lm≈0.85, rm≈0.21, midH=0.7 → panelW_conc ≈ 0.62*(6.5-1.76) ≈ 2.94
% panelH_conc = 2.94/4 ≈ 0.74  → figH ≈ bm + 2*0.74 + midV + tm ≈ 2.6
fig = journal_figure([6.5, 2.5], 2, 'tex');

axPi2  = axes(fig); hold(axPi2, 'on');
axTh0  = axes(fig); hold(axTh0, 'on');
axSpeed = axes(fig); hold(axSpeed, 'on'); box(axSpeed, 'on');

% ==========================================================================
% Panels (a) & (b): concentration field rendering
% ==========================================================================
configs = {axPi2, fullfile(baseDir,'center_speed_snapshot_thpi2.dat'), pi/2, c2, '{\it \theta} = {\it \pi}/2'; ...
           axTh0, fullfile(baseDir,'center_speed_snapshot_th0.dat'),   0,    c0, '{\it \theta} = 0'};

for j = 1:2
    ax    = configs{j,1};
    fname = configs{j,2};
    theta = configs{j,3};
    col   = configs{j,4};
    lbl   = configs{j,5};

    [xplt, yplt, phi, U, V] = load_snapshot(fname);

    imagesc(ax, xplt, yplt, phi);
    ax.YDir = 'normal';
    colormap(ax, cmocean('tempo'));
    clim(ax, [0 1]);
    axis(ax, 'equal'); axis(ax, 'tight');
    ax.XTick = []; ax.YTick = [];
    box(ax, 'on'); ax.XColor = 'k'; ax.YColor = 'k';

    % Circle overlays
    ct = cos(theta); st = sin(theta);
    for k = 1:Np
        xe = r*cos(phi_c) + xpos(k);
        ye = r*sin(phi_c);
        patch(ax, xe, ye, 'w', 'FaceColor', 'w', 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 2);
    end

    % Streamlines
    if j == 1
        % Panel (a): seeds equispaced on x-axis (y=0), forward + backward
        Nsl = 51;
        sx  = linspace(xplt(1), xplt(end), Nsl);
        sy  = zeros(1, Nsl);
        [X, Y] = meshgrid(xplt, yplt);
        h1 = streamline(ax, X, Y,  U,  V, sx, sy, [1 100]);
        h2 = streamline(ax, X, Y, -U, -V, sx, sy, [1 100]);
        set([h1; h2], 'Color', 'w', 'LineWidth', 1);
    else
        % Panel (b): seeds equispaced in y at each particle center x
        [X, Y] = meshgrid(xplt, yplt);

        % Seeds at each particle x, y in [-d/2, d/2]
        Nsl  = 3;
        sy_p = linspace(-d/3, 0, Nsl);
        sy_p = [sy_p, -sy_p];  % add seeds in upper strip
        Nsl = numel(sy_p);
        sx_p = repmat(xpos, 1, Nsl);
        sy_p = repelem(sy_p, 1, Np);

        % Seeds at center x, y in outer strips [ymin,-d/2] and [d/2,ymax]
        Nsl_out = 13; sc = 0.8;
        sy_out  = [linspace(yplt(1), -d*sc, Nsl_out), linspace(d*sc, yplt(end), Nsl_out)];
        sx_out  = zeros(1, numel(sy_out));
        sx_out = ones(1, numel(sy_out)) * xplt(1);  % all seeds at center x

        sx = [sx_p, sx_out];
        sy = [sy_p, sy_out];

        h1 = streamline(ax, X, Y,  U,  V, sx, sy, [1 1000]);
        h2 = streamline(ax, X, Y, -U, -V, sx, sy, [1 1000]);
        set([h1; h2], 'Color', 'w', 'LineWidth', 1);
    end
    
end

% Shared colorbar on axPi2
cb = colorbar(axPi2, 'Location', 'westoutside', 'Ticks', 0:0.5:1);

% ==========================================================================
% Panel (c): center speed vs N
% ==========================================================================
plot(axSpeed, N, uN0, 'o-', 'Color', c0, 'DisplayName', '{\it \theta} = 0');
plot(axSpeed, N, 1+0*N, '--', 'Color', c0, 'DisplayName', 'multipole');
plot(axSpeed, N, uN2, 'o-', 'Color', c2, 'DisplayName', '{\it \theta} = {\it \pi}/2');
plot(axSpeed, N, u_mp, '--', 'Color', c2, 'DisplayName', 'multipole');
% yline(axSpeed, 1, ':', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off', 'LineWidth', 2);

xlabel(axSpeed, 'number of particles, {\it N}');
ylabel(axSpeed, 'center speed, {\it u}_{cm} / {\it u}^{\it S}');
axSpeed.XTick = N';
legend(axSpeed, 'Location', 'northwest', 'Color', 'none');
xlim([0 14])
ylim(axSpeed, [0.8 1.4]);

% ==========================================================================
% Layout
% ==========================================================================
for ax = [axPi2, axTh0, axSpeed]
    ax.ActivePositionProperty = 'position';
    ax.Units = 'inches';
end
drawnow;

tiS = axSpeed.TightInset;
lm     = 0.8;
rm     = max(0.15, tiS(3)) + 0.06;
bm     = max(0.42, tiS(2)) + 0.08;
tm     = max(0.15, tiS(4)) + 0.06;
midH   = 1;    % horizontal gap between concentration strips and speed panel
midV   = 0.25;   % vertical gap between the two concentration strips

figW = fig.PaperSize(1);
figH = fig.PaperSize(2);

panelH_full = figH - bm - tm;

% Speed panel: right portion of figure
panelW_speed = 0.35 * (figW - lm - rm - midH);
panelW_conc  = figW - lm - rm - midH - panelW_speed;

% 4:1 aspect ratio for each concentration strip (matches the data grid)
panelH_conc = panelW_conc / 4;

xConc  = lm;
xSpeed = lm + panelW_conc + midH;

axPi2.Position  = [xConc,  bm + panelH_conc + midV, panelW_conc,  panelH_conc];
axTh0.Position  = [xConc,  bm,                       panelW_conc,  panelH_conc];
axSpeed.Position = [xSpeed, bm,                       panelW_speed, panelH_full];

drawnow;
axPi2.Position  = [xConc,  bm + panelH_conc + midV, panelW_conc,  panelH_conc];
axTh0.Position  = [xConc,  bm,                       panelW_conc,  panelH_conc];
axSpeed.Position = [xSpeed, bm,                       panelW_speed, panelH_full];

% Colorbar spanning both concentration panels
cb.Units = 'inches';
cbW = 0.1;   % colorbar width
cbG = 0.08;  % gap between colorbar and panels
cb.Position = [xConc - cbG - cbW, bm, cbW, 2*panelH_conc + midV];

cb.Title.String   = '{\it \phi}';
cb.Title.FontName = 'Times';
cb.Title.FontSize = axPi2.FontSize;

panel_label('\it (a)', axPi2);
panel_label('\it (b)', axTh0);
panel_label('\it (c)', axSpeed);

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
