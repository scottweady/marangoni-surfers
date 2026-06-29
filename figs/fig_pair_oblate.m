% Three-panel figure for ellipse pair (a=0.2, b=4.0).
% Left two: x-y trajectories for d=10 and d=20 with ellipse outlines.
% Right: separation distance vs time for each initial separation d.

close all;

baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir, '..', 'data', 'oblate');

if ~exist('P', 'var'), P = palette(); end

% ---- Load all separation data -----------------------------------------------
dfiles = dir(fullfile(dataDir, 'a0.2500_b4.0000_d*.dat'));
d_vals = zeros(1, numel(dfiles));
for k = 1:numel(dfiles)
    tok = sscanf(dfiles(k).name, 'a%f_b%f_d%f.dat');
    d_vals(k) = tok(3);
end
[d_vals, isort] = sort(d_vals);
dfiles = dfiles(isort);

% ---- Figure -----------------------------------------------------------------
fig = journal_figure([6.75, 2.25], 2, 'tex');

col1 = P.c1;
col2 = P.c2;

% --- Trajectory panels (d=10 and d=20) ---------------------------------------
traj_d = [12, 18];
axTraj = gobjects(1, 2);
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

    ax = subplot(1, 3, ip, 'Parent', fig);
    axTraj(ip) = ax;
    hold(ax, 'on'); box(ax, 'on'); axis(ax, 'equal');

    plot(ax, x1, y1, '-', 'Color', col1, 'DisplayName', 'particle 1');
    plot(ax, x2, y2, '-', 'Color', col2, 'DisplayName', 'particle 2');

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
        th = parts{p,3}(1); col = parts{p,4};
        ct = cos(th); st = sin(th);
        xe = ct*ex0 - st*ey0 + xc;
        ye = st*ex0 + ct*ey0 + yc;
        patch(ax, xe, ye, [0.92 0.92 0.92], 'EdgeColor', [0.2 0.2 0.2], ...
              'LineWidth', 1.5, 'HandleVisibility', 'off');
        xe_arc = ct*xa_act - st*ya_act + xc;
        ye_arc = st*xa_act + ct*ya_act + yc;
        plot(ax, xe_arc, ye_arc, '-', 'Color', col, 'LineWidth', 4, ...
             'HandleVisibility', 'off');
        quiver(ax, xc, yc, cos(th), sin(th), 3, 'Color', 'k', ...
               'LineWidth', 2, 'MaxHeadSize', 1, 'HandleVisibility', 'off');
    end

    xlabel(ax, '{\it x}');
    ylabel(ax, '{\it y}');
    axis(ax, 18 * [-1 1 -1 1]);
    xticks(ax, -20:10:20);
    yticks(ax, -20:10:20);
    % title(ax, sprintf('$d_0 = %g$', traj_d(ip)), 'Interpreter', 'latex');
end

set(axTraj(2), 'YTickLabel', [], 'YLabel', []);
% --- Right panel: separation distance vs time --------------------------------
ax2 = subplot(1, 3, 3, 'Parent', fig);
hold(ax2, 'on'); box(ax2, 'on');

cmap = cmocean('tempo', numel(dfiles) + 2);
cmap = cmap(2:end-1, :);
for k = 1:numel(dfiles)
    Tk = readmatrix(fullfile(dataDir, dfiles(k).name));
    tk  = Tk(:,1);
    sep = sqrt((Tk(:,5) - Tk(:,2)).^2 + (Tk(:,6) - Tk(:,3)).^2);
    plot(ax2, tk, sep, '-', 'Color', cmap(k,:), ...
         'DisplayName', strcat('{\it d}_0 = {}', num2str(d_vals(k))));
    if k == 1
        sep_inf = sep(end);
    end

end

plot(ax2, [0 5e4], [1 1] * sep_inf, ':', 'Color', 0.2 * [1 1 1], 'HandleVisibility', 'off');

xlabel(ax2, 'time, {\it t}');
ylabel(ax2, 'distance, {\it d}');
legend(ax2, 'Location', 'northeast', 'Color', 'none');

xlim([0 50000]);
ylim([8 32]);
yticks(8:8:32);

% --- Resize panels to fill figure (all in inches) ----------------------------
ml = 0.75;  mr = 0.20;  mb = 0.75;  mt = 0.30;  gap = 1.0;

fig.Units = 'inches';
figW = fig.Position(3);
figH = fig.Position(4);

sp1 = axTraj(1);
sp2 = axTraj(2);
sp3 = ax2;

set([sp1, sp2, sp3], 'Units', 'inches');

sp1.Position(1) = ml;
sp1.Position(2) = mb;
sp1.Position(3) = (figW - ml - mr - gap) / 3.25;
sp1.Position(4) = figH - mb - mt;

sp2.Position(1) = sp1.Position(1) + sp1.Position(3);
sp2.Position(2:4) = sp1.Position(2:4);

sp3.Position(1) = sp2.Position(1) + sp2.Position(3) + gap;
sp3.Position(3) = figW - sp3.Position(1) - mr;
sp3.Position([2 4]) = sp1.Position([2 4]);

% --- Panel labels ------------------------------------------------------------
panel_label('\it (a) d_{\rm0}\rm = 12', sp1);
panel_label('\it (b) d_{\rm0}\rm = 18', sp2);
panel_label('\it (c)', sp3);

function a = panel_label(str, ax)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', 20, 'Units', 'inches');
    a.Position(1) = ax.Position(1);
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4);
end
