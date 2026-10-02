% fig_expansion.m — render particle configurations from data/frame_XXXXXX.dat

close all

a = sqrt(1.0);
b = 1.0 / a;
U = 0.098175;

% ── Load frames ───────────────────────────────────────────────────────────────

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'patch');
files   = dir(fullfile(dataDir, 'frame_*.dat'));
files   = sort({files.name});
Nframes = numel(files);

frames = cell(Nframes, 1);
for k = 1:Nframes
    frames{k} = load(fullfile(dataDir, files{k}));  % columns: t  x_re  x_im  θ  vx  vy  ω
end

% ── Fixed axes from final frame ───────────────────────────────────────────────

last  = frames{end};
pad   = 2 * max(a, b);
xlims = [min(last(:,2))-pad,  max(last(:,2))+pad];
ylims = [min(last(:,3))-pad,  max(last(:,3))+pad];

% ── Ellipse template (local frame) ───────────────────────────────────────────

s  = linspace(0, 2*pi, 64);
ex = a * cos(s);
ey = b * sin(s);

% ── Precompute radius and speeds for all frames ───────────────────────────────

Np = size(frames{1}, 1);
tv = zeros(Nframes, 1);
Rv = zeros(Nframes, 1);
speeds = cell(Nframes, 1);   % speed of each particle at each frame

for k = 1:Nframes
    d  = frames{k};
    xbar = mean(d(:,2));
    ybar = mean(d(:,3));
    tv(k) = d(1,1);
    Rv(k) = max(hypot(d(:,2) - xbar, d(:,3) - ybar));
end

for k = 1:Nframes
    d  = frames{k};
    xbar = mean(d(:,2));
    ybar = mean(d(:,3));
    rx = d(:,2) - xbar;
    ry = d(:,3) - ybar;
    r  = hypot(rx, ry);
    speeds{k} = (d(:,5).*rx + d(:,6).*ry) ./ r;   % radial component of velocity
end

Ubar = cellfun(@mean, speeds);

% ── Figure ───────────────────────────────────────────────────────────────────

fig = journal_figure([6.5 2.5], 2, 'tex');

tf = tv(end);

if ~exist('P', 'var'), P = palette(); end
dist_cmap = P.snap;
line_col  = P.col;

% Col 1a: radius vs time
sp_R = axes(fig); hold(sp_R, 'on'); box(sp_R, 'on');
plot(sp_R, tv, Rv/Rv(1), '-', 'Color', line_col);
ylabel(sp_R, '{\it R}/{\it R}_0');
xlim(sp_R, [0, tf]);
ylim(sp_R, [0, max(Rv)/Rv(1)*1.1]);
set(sp_R, 'XTickLabel', {});

% Col 1b: mean speed vs time
sp_V = axes(fig); hold(sp_V, 'on'); box(sp_V, 'on');
plot(sp_V, tv, Ubar/U, '-', 'Color', line_col);
xlabel(sp_V, 'time, {\it t}');
ylabel(sp_V, '\langle |{\it\bf u}_{cm}| \rangle / |{\it\bf u}_{cm,0}|');
xlim(sp_V, [0, tf]);
ylim(sp_V, [1 1.5]);
yticks([1 1.25 1.5]);

% Col 2: particle configurations at t=0 and t=final
sp_C = axes(fig); hold(sp_C, 'on'); axis(sp_C, 'equal');
box(sp_C, 'on'); set(sp_C, 'XTick', [], 'YTick', []);
d_ring = 5.0;   % inter-ring spacing
snap_cols = [dist_cmap(1,:); dist_cmap(3,:)];

% Shaded extents (drawn first so particles appear on top)
th_ring  = linspace(0, 2*pi, 256);
d0_      = frames{1};       dN_      = frames{end};
xbar0_     = mean(d0_(:,2));  ybar0_     = mean(d0_(:,3));
xbarN_     = mean(dN_(:,2));  ybarN_     = mean(dN_(:,3));
R0_max   = max(hypot(d0_(:,2)-xbar0_, d0_(:,3)-ybar0_));
RN_min   = min(hypot(dN_(:,2)-xbarN_, dN_(:,3)-ybarN_));
RN_max   = max(hypot(dN_(:,2)-xbarN_, dN_(:,3)-ybarN_));

fill(sp_C, xbar0_ + R0_max*cos(th_ring), ybar0_ + R0_max*sin(th_ring), ...
     snap_cols(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);

for sgn = [1, -1]
    th_half = linspace(0, sgn*pi, 128);
    xo_ = xbarN_ + RN_max*cos(th_half);
    yo_ = ybarN_ + RN_max*sin(th_half);
    xi_ = xbarN_ + RN_min*cos(fliplr(th_half));
    yi_ = ybarN_ + RN_min*sin(fliplr(th_half));
    fill(sp_C, [xo_, xi_], [yo_, yi_], snap_cols(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
end

for ki = [1, Nframes]
    fdata = frames{ki};
    x   = fdata(:,2);  y = fdata(:,3);  th = fdata(:,4);
    col = snap_cols(1 + (ki > 1), :);
    for i = 1:numel(x)
        R  = [cos(th(i)) -sin(th(i)); sin(th(i)) cos(th(i))];
        xy = R * [ex; ey];
        fill(sp_C, x(i)+xy(1,:), y(i)+xy(2,:), col, 'EdgeColor', 'none');
    end
end
xpad = 0.05 * diff(xlims);  ypad = 0.05 * diff(ylims);
xlim(sp_C, xlims + [-xpad, xpad]); ylim(sp_C, ylims + [-ypad, ypad]);

% Scale bar: 50 particle radii
sb_len = 25 * a;
sb_x0  = xlims(1) + 0*diff(xlims);
sb_y   = ylims(1) + 0*diff(ylims);
line(sp_C, [sb_x0, sb_x0+sb_len], [sb_y, sb_y], 'Color', 'k', 'LineWidth', 2);
text(sp_C, sb_x0 + sb_len/2, sb_y + 0.02*diff(ylims), '25{\it R}', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'FontSize', sp_C.FontSize);

% Arrows at exact angles 0, pi/4, ..., 7pi/4
d0  = frames{1};    dN = frames{end};
xbar0 = mean(d0(:,2)); ybar0 = mean(d0(:,3));
xbarN = mean(dN(:,2)); ybarN = mean(dN(:,3));
R1  = max(hypot(d0(:,2)-xbar0, d0(:,3)-ybar0));
R2  = min(hypot(dN(:,2)-xbarN, dN(:,3)-ybarN));
phi = (0:7)' * pi/4;
dR = 5;
x0v = xbar0 + (R1+dR)*cos(phi);  y0v = ybar0 + (R1+dR)*sin(phi);
xNv = xbarN + (R2-dR)*cos(phi);  yNv = ybarN + (R2-dR)*sin(phi);
h = quiver(sp_C, x0v, y0v, xNv-x0v, yNv-y0v, 0, 'Color', [0.5 0.5 0.5]);
h.Color = 'k';  h.LineWidth = 1;  h.MaxHeadSize = 0.1;

% Col 3: surfactant concentration snapshot (from expansion_snapshot.jl)
fid  = fopen(fullfile(dataDir, 'expansion_snapshot.dat'));
hdr  = fgetl(fid);
Nx   = str2double(regexp(hdr, 'Nx=(\d+)', 'tokens', 'once'));
Ny   = str2double(regexp(hdr, 'Ny=(\d+)', 'tokens', 'once'));
t_S  = str2double(regexp(hdr, 't=([\d.]+)', 'tokens', 'once'));
xplt = str2num(fgetl(fid));  %#ok<ST2NM>
yplt = str2num(fgetl(fid));  %#ok<ST2NM>
phiS = fscanf(fid, '%f', [Nx, Ny])';
fclose(fid);

sp_S = axes(fig); hold(sp_S, 'on'); box(sp_S, 'on');
imagesc(sp_S, xplt, yplt, phiS, 'AlphaData', ~isnan(phiS));
sp_S.YDir = 'normal';
axis(sp_S, 'equal'); axis(sp_S, 'tight');
set(sp_S, 'XTick', [], 'YTick', []);
colormap(sp_S, P.cmap());
clim(sp_S, [0 1]);
cb = colorbar(sp_S, 'Location', 'eastoutside', 'Ticks', 0:0.5:1);
cb.Title.String = '\phi';

[~, kS] = min(abs(tv - t_S));
dS = frames{kS};
for i = 1:size(dS, 1)
    R  = [cos(dS(i,4)) -sin(dS(i,4)); sin(dS(i,4)) cos(dS(i,4))];
    xy = R * [ex; ey];
    fill(sp_S, dS(i,2)+xy(1,:), dS(i,3)+xy(2,:), 'w', 'EdgeColor', 'k', 'LineWidth', 0.5);
end
text(sp_S, 0.97, 0.03, sprintf('{\\it t} = %.0f', t_S), 'Units', 'normalized', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
     'Color', 'w', 'FontSize', sp_S.FontSize);

% ── Layout ────────────────────────────────────────────────────────────────────

lm = 0.2;   % left margin
rm = 0.80;   % right margin (room for colorbar)
bm = 0.80;   % bottom margin
tm = 0.0;   % top margin
cg = 1.00;   % gap between columns
vg = 0.2;   % vertical gap between stacked left panels

figW = fig.PaperSize(1);
figH = fig.PaperSize(2);

panelW  = (figW - lm - rm - 2*cg) / 3;
panelH  = figH - bm - tm;

% Height of (a) after axis equal: data fills panelW in x, height follows
height_C = diff(ylims) / diff(xlims) * panelW;
halfH    = (height_C - vg) / 2;

col1x = lm;
col2x = lm +   panelW + cg;
col3x = lm + 2*panelW + 2*cg;

for ax = [sp_R, sp_V, sp_C, sp_S]
    ax.Units = 'inches';
    ax.ActivePositionProperty = 'position';
end

sp_C.Position = [col1x, bm,              panelW, height_C];
sp_R.Position = [col2x, bm + halfH + vg, panelW, halfH];
sp_V.Position = [col2x, bm,              panelW, halfH];

sp_S.Position = [col3x, bm,              panelW, height_C];

% ── Panel labels ──────────────────────────────────────────────────────────────

panel_label('(a)', sp_C);
panel_label('(b)', sp_R);
panel_label('(c)', sp_V);
panel_label('(d)', sp_S);

function panel_label(str, ax)
    text(ax, 0.01, 0.99, ['\it ' str], 'Units', 'normalized', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
         'FontSize', ax.FontSize);
end
