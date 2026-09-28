% fig_suspension_config.m — five snapshots of the particle configuration
% over time for three ellipse aspect ratios (scripts/suspension/suspension.jl
% output), one snapshot color-scale per panel.

close all;

baseDir   = fullfile(fileparts(mfilename('fullpath')), '..', 'scripts', 'suspension');
ars       = [0.5, 1.0, 2.0];
Nsnap     = 5;
alpha_arc = pi / 2;   % active (Neumann) arc half-angle, matches suspension.jl
P         = palette();
cmap      = cmocean('tempo', 5);

data = cell(1, numel(ars));
for k = 1:numel(ars)
    frameDir = fullfile(baseDir, sprintf('frames_%.1f', ars(k)));
    files    = dir(fullfile(frameDir, 'frame_*.dat'));
    names    = sort({files.name});

    t = zeros(1, numel(names));
    x = cell(1, numel(names));
    y = cell(1, numel(names));
    yspan = zeros(1, numel(names));
    for j = 1:numel(names)
        A    = readmatrix(fullfile(frameDir, names{j}));
        t(j) = A(1, 1);
        x{j} = A(:, 2);
        y{j} = A(:, 3);
        yspan(j) = max(y{j}) - min(y{j});
    end

    stop = find(yspan >= 12, 1, 'first');
    if isempty(stop), stop = numel(t); end
    t = t(1:stop); x = x(1:stop); y = y(1:stop);

    idx     = round(linspace(1, numel(t), Nsnap));
    data{k} = {t(idx), x(idx), y(idx)};
end

fig = journal_figure([6.75, 2], 2, 'tex');
ax  = gobjects(1, numel(ars));
cb  = gobjects(1, numel(ars));

for k = 1:numel(ars)
    ax(k) = axes(fig); hold(ax(k), 'on'); box(ax(k), 'on'); axis(ax(k), 'equal');

    [t, x, y] = data{k}{:};
    tmax = max(t);
    for j = 1:numel(t)
        c    = cmap(j,:);
        yoff = y{j} - min(y{j});
        plot(ax(k), x{j}, yoff, '-o', 'Color', c, 'MarkerSize', 5, ...
            'MarkerFaceColor', 'none');
    end

    xlabel(ax(k), '{\it x}');
    xlim(45 * [-1 1])
    ylim(40 * [-0.2 1])

    schem_scale = 7;
    a = schem_scale * sqrt(ars(k)); b = schem_scale / sqrt(ars(k));
    xl = xlim(ax(k)); yl = ylim(ax(k));
    xc = xl(2) - 0.10 * diff(xl) - 4;
    yc = yl(2) - 0.10 * diff(yl) - 7;
    draw_particle_schematic(ax(k), xc, yc, a, b, pi/2, alpha_arc, P.col);

    colormap(ax(k), cmap);
    cb(k) = colorbar(ax(k), 'northoutside');
    clim(ax(k), [0 tmax]);
    cb(k).Title.String = '{\it t}';

    if k > 1
        set(ax(k), 'YTickLabel', {});
    end
end
ylabel(ax(1), '{\it y}');
linkaxes(ax, 'y');

% ── Layout ────────────────────────────────────────────────────────────────────
fig.Units = 'inches'; figW = fig.Position(3); figH = fig.Position(4);

ax(1).Units = 'inches';
ax(2).Units = 'inches';
ax(3).Units = 'inches';

ax(1).Position(1) = 0;
ax(1).Position(2) = 0.8;
ax(1).Position(4) = 2.1;

ax(2).Position(1) = sum(ax(1).Position([1 3])) - 1;
ax(2).Position(2:4) = ax(1).Position(2:4);

ax(3).Position(1) = sum(ax(2).Position([1 3])) - 1;
ax(3).Position(2:4) = ax(1).Position(2:4);

for i = 1:3
    cb(i).Units = 'inches';
    cb(i).Position(2) = sum(ax(i).Position([2 4])) + 0.1;
    cb(i).Position(3) = 0.75 * ax(i).Position(3);
    cb(i).Position(1) = ax(i).Position(1) + 0.5 * ax(i).Position(3) - 0.5 * cb(i).Position(3);
    cb(i).Position(4) = 0.05 * cb(i).Position(3);
end

% cbY = ax(1).Position(2) + ax(1).Position(4) + 0.15;
% cbH = 0.2;
% for k = 1:numel(ars)
%     cb(k).Units = 'inches';
%     cb(k).Position = [ax(k).Position(1), cbY, ax(k).Position(3), cbH];
% end

labels = {'\it (a)', '\it (b)', '\it (c)'};
for k = 1:numel(ars)
    panel_label(labels{k}, ax(k), 0.55, -0.01);
end

function a = panel_label(str, ax, xOffset, yOffset)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', ax.FontSize, 'Units', 'inches');
    a.Position(1) = ax.Position(1) + xOffset;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + yOffset;
end

function draw_particle_schematic(ax, xc, yc, a, b, theta, alpha_arc, col)
    t     = linspace(0, 2 * pi, 200);
    t_act = linspace(alpha_arc, 2 * pi - alpha_arc, 200);
    ct = cos(theta); st = sin(theta);

    ex = a * cos(t) * ct - b * sin(t) * st + xc;
    ey = a * cos(t) * st + b * sin(t) * ct + yc;
    fill(ax, ex, ey, 0.9 * [1 1 1], 'EdgeColor', 'k', 'LineWidth', 2);

    ex_act = a * cos(t_act) * ct - b * sin(t_act) * st + xc;
    ey_act = a * cos(t_act) * st + b * sin(t_act) * ct + yc;
    plot(ax, ex_act, ey_act, '-', 'Color', col, 'LineWidth', 4);
end
