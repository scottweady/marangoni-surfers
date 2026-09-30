% fig_pair_antiparallel.m — anti-parallel ellipse pairs: (a)-(c) concentration + velocity
% fields at d0 = 8 across aspect ratio; (d),(e) x-y trajectories of an oblate pair
% (a=0.25, b=4) for d0 = 12, 18; (f) separation distance vs time for each d0.

close all;

baseDir = fileparts(mfilename('fullpath'));
P = palette();

fig = journal_figure([6.75, 4.35], 2, 'tex');

% ── Fields ────────────────────────────────────────────────────────────────────
dataDir = fullfile(baseDir, '..', 'data', 'antiparallel');
col   = P.col;
aspan = [1.0/2.0, 1.0, 2.0];
th    = [-pi/2, pi/2];
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
    colormap(ax(k), P.cmap());
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

% ── Oblate pair: load all separation data ─────────────────────────────────────
dataDir = fullfile(baseDir, '..', 'data', 'oblate');

dfiles = dir(fullfile(dataDir, 'a0.2500_b4.0000_d*.dat'));
d_vals = zeros(1, numel(dfiles));
for k = 1:numel(dfiles)
    tok = sscanf(dfiles(k).name, 'a%f_b%f_d%f.dat');
    d_vals(k) = tok(3);
end
[d_vals, isort] = sort(d_vals);
dfiles = dfiles(isort);

col1 = P.c1;
col2 = P.c2;

% --- Trajectory panels -------------------------------------------------------
traj_d = [12, 18];
axTraj = [axes(fig), axes(fig)];
for ip = 1:2
    fname_traj = sprintf('a0.2500_b4.0000_d%g.0.dat', traj_d(ip));
    T  = readmatrix(fullfile(dataDir, fname_traj));
    id = 1 : min(400, size(T,1));
    T  = T(id, :);
    x1 = T(:,2);  y1 = T(:,3);  th1 = T(:,4);
    x2 = T(:,5);  y2 = T(:,6);  th2 = T(:,7);

    tokens = sscanf(fname_traj, 'a%f_b%f_d%f.dat');
    a_ell  = tokens(1);
    b_ell  = tokens(2);

    axk = axTraj(ip);
    hold(axk, 'on'); box(axk, 'on'); axis(axk, 'equal');

    plot(axk, x1, y1, '-', 'Color', col1, 'DisplayName', 'particle 1');
    plot(axk, x2, y2, '-', 'Color', col2, 'DisplayName', 'particle 2');

    phi_e  = linspace(0, 2*pi, 200);
    ex0    = a_ell * cos(phi_e);
    ey0    = b_ell * sin(phi_e);
    alp    = pi/2;
    phi_act = linspace(alp, 2*pi - alp, 150);
    xa_act  = a_ell * cos(phi_act);
    ya_act  = b_ell * sin(phi_act);

    parts = {x1, y1, th1, col1; x2, y2, th2, col2};
    for p = 1:2
        xc = parts{p,1}(1); yc = parts{p,2}(1);
        thp = parts{p,3}(1); colp = parts{p,4};
        ct = cos(thp); st = sin(thp);
        xe = ct*ex0 - st*ey0 + xc;
        ye = st*ex0 + ct*ey0 + yc;
        patch(axk, xe, ye, [0.92 0.92 0.92], 'EdgeColor', [0.2 0.2 0.2], ...
              'LineWidth', 1.5, 'HandleVisibility', 'off');
        xe_arc = ct*xa_act - st*ya_act + xc;
        ye_arc = st*xa_act + ct*ya_act + yc;
        plot(axk, xe_arc, ye_arc, '-', 'Color', colp, 'LineWidth', 4, ...
             'HandleVisibility', 'off');
        quiver(axk, xc, yc, cos(thp), sin(thp), 3, 'Color', 'k', ...
               'LineWidth', 2, 'MaxHeadSize', 1, 'HandleVisibility', 'off');
    end

    xlabel(axk, '{\it x}');
    ylabel(axk, '{\it y}');
    axis(axk, 18 * [-1 1 -1 1]);
    xticks(axk, -20:10:20);
    yticks(axk, -20:10:20);
end

set(axTraj(2), 'YTickLabel', [], 'YLabel', []);

% --- Separation distance vs time ---------------------------------------------
axDist = axes(fig);
hold(axDist, 'on'); box(axDist, 'on');

cmap = P.cmap(numel(dfiles) + 2);
cmap = cmap(2:end-1, :);
for k = 1:numel(dfiles)
    Tk = readmatrix(fullfile(dataDir, dfiles(k).name));
    tk  = Tk(:,1);
    sep = sqrt((Tk(:,5) - Tk(:,2)).^2 + (Tk(:,6) - Tk(:,3)).^2);
    plot(axDist, tk, sep, '-', 'Color', cmap(k,:), ...
         'DisplayName', strcat('{\it d}_0 = {}', num2str(d_vals(k))));
    if k == 1
        sep_inf = sep(end);
    end
end

plot(axDist, [0 5e4], [1 1] * sep_inf, ':', 'Color', 0.2 * [1 1 1], 'HandleVisibility', 'off');

xlabel(axDist, 'time, {\it t}');
ylabel(axDist, 'distance, {\it d}');
legend(axDist, 'Location', 'northeast', 'Color', 'none');

xlim(axDist, [0 50000]);
ylim(axDist, [8 32]);
yticks(axDist, 8:8:32);

% ── Layout ────────────────────────────────────────────────────────────────────
topH = 4.2; botH = 4.5;

fig.Units = 'inches'; figW = fig.Position(3);
set([ax, axTraj, axDist], 'Units', 'inches');
set(ax, 'ActivePositionProperty', 'position');

% fields row
ml = 0.75; mr = 0.75; mb = 0.25; mt = 0.15; gap = 0.1;
panelW = (figW - ml - mr - (Np-1)*gap) / Np;
panelH = topH - mb - mt;
for k = 1:Np
    ax(k).Position = [ml + (k-1)*(panelW+gap), botH + mb, panelW, panelH];
end

% oblate row
ml = 0.75; mr = 0.75; mb = 0.75; mt = 0.30; gap = 1.0;
trajW  = (figW - ml - mr - gap) / 3.25;
panelH = botH - mb - mt;
axTraj(1).Position = [ml,         mb, trajW, panelH];
axTraj(2).Position = [ml + trajW, mb, trajW, panelH];
axDist.Position    = [ml + 2*trajW + gap, mb, figW - ml - 2*trajW - gap - mr, panelH];

% ── Panel labels ──────────────────────────────────────────────────────────────
letters = 'abc';
for k = 1:Np
    a = aspan(k); b = 1.0/a;
    panel_label(sprintf('\\it (%s) {\\it a/b} \\rm = %.2f', letters(k), a/b), ax(k), -0.2, -0.1);
end
panel_label('\it (d) d_{\rm0}\rm = 12', axTraj(1), 0, 0);
panel_label('\it (e) d_{\rm0}\rm = 18', axTraj(2), 0, 0);
panel_label('\it (f)', axDist, 0, 0);

function a = panel_label(str, ax, xOffset, yOffset)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', 20, 'Units', 'inches');
    a.Position(1) = ax.Position(1) + xOffset;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + yOffset;
end
