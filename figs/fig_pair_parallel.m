% fig_pair_parallel.m — parallel ellipse pairs: (a)-(c) concentration + velocity
% fields at d0 = 8, (d) separation and (e) orientation vs time across aspect ratio.

close all;

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'parallel');
P    = palette();
col  = P.col;

fig = journal_figure([6.75, 4.35], 2, 'tex');

% ── Fields ────────────────────────────────────────────────────────────────────
aspan = [1.0/2.0, 1.0, 2.0];
th    = [pi/2, pi/2];
Delta = 8.0;
alpha_arc = pi/2;
Np    = numel(aspan);

ax = arrayfun(@(~) axes(fig), 1:Np);

for k = 1:Np
    a = aspan(k);  b = 1.0/a;
    ax(k).FontSize = 20;
    fname = sprintf('fields_a%.4f_b%.4f.dat', a, b);
    fid   = fopen(fullfile(dataDir, fname));
    hdr   = fgetl(fid);
    Nx    = str2double(regexp(hdr, 'Nx=(\d+)', 'tokens', 'once'));
    Ny    = str2double(regexp(hdr, 'Ny=(\d+)', 'tokens', 'once'));
    xplt  = str2num(fgetl(fid));  %#ok<ST2NM>
    yplt  = str2num(fgetl(fid));  %#ok<ST2NM>
    phi   = fscanf(fid, '%f', [Nx, Ny])';
    U     = fscanf(fid, '%f', [Nx, Ny])';
    V     = fscanf(fid, '%f', [Nx, Ny])';
    fclose(fid);

    hold(ax(k), 'on'); box(ax(k), 'on');
    set(ax(k), 'BoxStyle', 'full', 'LineWidth', 1.5);
    imagesc(ax(k), xplt, yplt, phi);
    ax(k).YDir = 'normal';
    axis(ax(k), 'equal'); axis(ax(k), 'tight');
    colormap(ax(k), cmocean('tempo'));
    clim(ax(k), [0 1]);

    ax(k).XLabel.String = '{\it x}';
    if k == 1
        ax(k).YLabel.String = '{\it y}';
    else
        ax(k).YTick = [];
    end

    xcm  = [-Delta/2, Delta/2];
    ycm  = [0, 0];
    phi_e   = linspace(0, 2*pi, 200);
    ex0     = a*cos(phi_e); ey0 = b*sin(phi_e);
    phi_act = linspace(alpha_arc, 2*pi - alpha_arc, 200);
    ex_act  = a*cos(phi_act); ey_act = b*sin(phi_act);
    for j = 1:2
        ct = cos(th(j)); st = sin(th(j));
        patch(ax(k), ct*ex0 - st*ey0 + xcm(j), st*ex0 + ct*ey0 + ycm(j), ...
              'w', 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5);
        plot(ax(k), ct*ex_act - st*ey_act + xcm(j), st*ex_act + ct*ey_act + ycm(j), ...
             '-', 'Color', col, 'LineWidth', 2);
    end

    before = ax(k).Children;
    streamslice(ax(k), xplt, yplt, U, V, 2);
    h = setdiff(ax(k).Children, before);
    set(h, 'Color', 'w', 'LineWidth', 1);
end

cb = colorbar(ax(end), 'Location', 'eastoutside', 'Ticks', 0:0.5:1);
cb.Title.String = '\phi';

% ── Separation and orientation vs time ───────────────────────────────────────
files = dir(fullfile(dataDir, 'a*_d*.dat'));
a_vals = zeros(1, numel(files));
b_vals = zeros(1, numel(files));
for k = 1:numel(files)
    tok = sscanf(files(k).name, 'a%f_b%f_d%f.dat');
    a_vals(k) = tok(1);
    b_vals(k) = tok(2);
end
[~, isort] = sort(a_vals);
files  = files(isort);
a_vals = a_vals(isort);
b_vals = b_vals(isort);

cmap = colorcet('D11', 'N', numel(files)).^1.5;

axSep = axes(fig); hold(axSep, 'on'); box(axSep, 'on');
axOmg = axes(fig); hold(axOmg, 'on'); box(axOmg, 'on');

for k = 1:numel(files)
    A    = readmatrix(fullfile(files(k).folder, files(k).name));
    A    = A(all(~isnan(A), 2), :);
    c    = cmap(k, :);
    lbl  = ['{\it a/b} = ' num2str(a_vals(k)/b_vals(k), '%.2f')];
    t    = A(:,1);
    sep  = A(:,8);
    ang1 = unwrap(A(:,4));
    ang2 = unwrap(A(:,7));
    plot(axSep, t, sep/sep(1), 'Color', c, 'DisplayName', lbl);
    plot(axOmg, t, ang1 - pi/2, 'Color', c, 'HandleVisibility', 'off');
    plot(axOmg, t, ang2 - pi/2, 'Color', c, 'LineStyle', '--', 'HandleVisibility', 'off');

    if a_vals(k) > b_vals(k)
        plot(axSep, t(end), sep(end)/sep(1), 'o', 'Color', c, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
        plot(axOmg, t(end), ang1(end) - pi/2, 'o', 'Color', c, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
        plot(axOmg, t(end), ang2(end) - pi/2, 'o', 'Color', c, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    end
end

linkaxes([axSep, axOmg], 'x');
xlabel(axSep, 'time, {\it t}');
xlabel(axOmg, 'time, {\it t}');
ylabel(axSep, 'dist({\it\Gamma}_1, {\it\Gamma}_2)');
ylabel(axOmg, '{\it \theta} - \pi/2');
legend(axSep, 'Location', 'southwest', 'Color', 'none');

xlim(axSep, [0 500]); xlim(axOmg, [0 500]);

% ── Layout ────────────────────────────────────────────────────────────────────
topH = 4.2; botH = 4.5;

fig.Units = 'inches'; figW = fig.Position(3);
set([ax, axSep, axOmg], 'Units', 'inches');
set(ax, 'ActivePositionProperty', 'position');

% fields row
ml = 0.75; mr = 0.75; mb = 0.25; mt = 0.15; gap = 0.1;
panelW = (figW - ml - mr - (Np-1)*gap) / Np;
panelH = topH - mb - mt;
for k = 1:Np
    ax(k).Position = [ml + (k-1)*(panelW+gap), botH + mb, panelW, panelH];
end

% time series row
ml = 1.0; mr = 0.75; mb = 0.75; mt = 0.5; gapH = 1.0;
panelW = (figW - ml - mr - gapH) / 2;
panelH = botH - mb - mt;
axSep.Position = [ml,                 mb, panelW, panelH];
axOmg.Position = [ml + panelW + gapH, mb, panelW, panelH];

letters = 'abc';
for k = 1:Np
    a = aspan(k); b = 1.0/a;
    panel_label(sprintf('\\it (%s) {\\it a/b} \\rm = %.2f', letters(k), a/b), ax(k), -0.2, -0.1);
end
panel_label('\it (d)', axSep, -0.95, 0.2);
panel_label('\it (e)', axOmg, -0.95, 0.2);

function a = panel_label(str, ax, xOffset, yOffset)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', ax.FontSize, 'Units', 'inches');
    a.Position(1) = ax.Position(1) + xOffset;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + yOffset;
end
