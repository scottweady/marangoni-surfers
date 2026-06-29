close all;

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'prolate');
P    = palette();
cmap = P.lines;
col  = P.col;

fig    = journal_figure([6.5, 2.75], 2, 'tex');
axConc = axes(fig); hold(axConc, 'on'); box(axConc, 'on'); axConc.XColor = 'k'; axConc.YColor = 'k';
axSep  = axes(fig); hold(axSep,  'on'); box(axSep,  'on');
axOmg  = axes(fig); hold(axOmg,  'on'); box(axOmg,  'on');

% --- Panel (a): concentration field ------------------------------------------
fid  = fopen(fullfile(dataDir, 'concentration_snapshot.dat'));
hdr  = fgetl(fid);
Nx   = str2double(regexp(hdr, 'Nx=(\d+)', 'tokens', 'once'));
Ny   = str2double(regexp(hdr, 'Ny=(\d+)', 'tokens', 'once'));
xplt = str2num(fgetl(fid));  %#ok<ST2NM>
yplt = str2num(fgetl(fid));  %#ok<ST2NM>
phi  = fscanf(fid, '%f', [Nx, Ny])';
U    = fscanf(fid, '%f', [Nx, Ny])';
V    = fscanf(fid, '%f', [Nx, Ny])';
fclose(fid);

imagesc(axConc, xplt, yplt, phi);
axConc.YDir = 'normal';
axis(axConc, 'equal'); axis(axConc, 'tight');
colormap(axConc, cmocean('tempo'));
axConc.XTick = []; axConc.YTick = [];

B   = readmatrix(fullfile(dataDir, 'a2.0000_b0.5000_d16.0.dat'));
row = B(end-15, :);
xcm_snap = [row(2), row(5)];
ycm_snap = [row(3), row(6)];
th_snap  = [row(4), row(7)];
phi_e = linspace(0, 2*pi, 200);
ex0 = 2.0*cos(phi_e); ey0 = 0.5*sin(phi_e);
phi_act = linspace(pi/2, 3*pi/2, 200);
ex_act = 2.0*cos(phi_act); ey_act = 0.5*sin(phi_act);
for k = 1:2
    ct = cos(th_snap(k)); st = sin(th_snap(k));
    patch(axConc, ct*ex0 - st*ey0 + xcm_snap(k), st*ex0 + ct*ey0 + ycm_snap(k), ...
          'w', 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 2);
    plot(axConc, ct*ex_act - st*ey_act + xcm_snap(k), st*ex_act + ct*ey_act + ycm_snap(k), ...
         '-', 'Color', col, 'LineWidth', 3);
end
xs_seed = linspace(xplt(1), xplt(end), 12);
ys_pos  = linspace(yplt(1), yplt(end), 4);
ys_pos = [ys_pos, 47 47 48 48];
[ysg,xsg] = meshgrid(ys_pos, xs_seed);
xy_fwd = stream2(xplt, yplt,  U,  V, xsg, ysg);
xy_bwd = stream2(xplt, yplt, -U, -V, xsg, ysg);
h = streamline(axConc, [xy_fwd; xy_bwd]);
set(h, 'Color', 'w', 'LineWidth', 1.5);

cb       = colorbar(axConc, 'Location', 'westoutside', 'Ticks', 0:0.5:1);
clim(axConc, [0 1]);

% --- Panels (b,c): separation and angular velocity vs time -------------------
files = dir(fullfile(dataDir, '*.dat'));
for k = 1:numel(files)
    a = parse_param(files(k).name, 'a');
    b = parse_param(files(k).name, 'b');
    if isnan(a); continue; end
    A   = readmatrix(fullfile(files(k).folder, files(k).name));
    A   = A(all(~isnan(A), 2), :);
    c   = cmap(round(a),:);
    lbl = ['{\it a/b} = ' num2str(round(a/b))];
    sep  = A(:,8);
    t    = A(:,1);
    ang1 = unwrap(A(:,4));
    ang2 = unwrap(A(:,7));
    plot(axSep, t, sep,  'Color', c, 'DisplayName', lbl);
    plot(axOmg, t, ang1-pi/2, 'Color', c, 'HandleVisibility', 'off');
    plot(axOmg, t, ang2-pi/2, 'Color', c, 'LineStyle', '--', 'HandleVisibility', 'off');
    plot(axSep, t(end), sep(end),  'o', 'Color', c, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    plot(axOmg, t(end), ang1(end)-pi/2, 'o', 'Color', c, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    plot(axOmg, t(end), ang2(end)-pi/2, 'o', 'Color', c, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
end

linkaxes([axSep, axOmg], 'x');
axSep.XTickLabel = {};
ylabel(axSep, 'dist({\it\Gamma}_1, {\it\Gamma}_2)');
ylabel(axOmg, '{\it \theta} - \pi/2');
xlabel(axOmg, 'time, {\it t}');
legend(axSep, 'Location', 'northeast', 'Color', 'none');

% --- Layout ------------------------------------------------------------------
ml = 1.0; mr = 0.20; mb = 0.75; mt = 0.30;
gapH = 1.25;   % horizontal gap between (a) and right column
gapV = 0.35;   % vertical gap between (b) and (c)

fig.Units = 'inches'; figW = fig.Position(3); figH = fig.Position(4);
set([axConc, axSep, axOmg], 'Units', 'inches');
panelW  = (figW - ml - mr - gapH) / 2;
panelH  = figH - mb - mt;
panelHR = (panelH - gapV) / 2;

xR = ml + panelW + gapH;
axConc.Position = [ml, mb,                   panelW, panelH];
axSep.Position  = [xR, mb + panelHR + gapV,  panelW, panelHR];
axOmg.Position  = [xR, mb,                   panelW, panelHR];

cb.Units = 'inches';
cb.Title.String = '\phi';


set(axSep, 'FontSize', axConc.FontSize);
set(axOmg, 'FontSize', axConc.FontSize);

panel_label('\it (a)', axConc, -1.0, 0.2);
panel_label('\it (b)', axSep,  -0.9, 0.2);
panel_label('\it (c)', axOmg,  -0.9, 0.2);

% --- Helpers -----------------------------------------------------------------
function val = parse_param(fname, key)
    tok = regexp(fname, ['(?:^|_)' key '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'], 'tokens', 'once');
    if isempty(tok); val = NaN; else; val = str2double(tok{1}); end
end

function a = panel_label(str, ax, xOffset, yOffset)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', ax.FontSize, 'Units', 'inches');
    a.Position(1) = ax.Position(1) + xOffset;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + yOffset;
end
