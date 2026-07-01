% fig_pair_conc.m — three-panel concentration snapshots: disk / prolate / oblate

close all;

P   = palette();
col = P.col;

fig = journal_figure([6.75, 2.25], 2, 'tex');

cases = struct( ...
    'dir',     {'disk',                        'prolate',                      'oblate'}, ...
    'a',       {1.0,                            2.0,                            0.25}, ...
    'b',       {1.0,                            0.5,                            4.0}, ...
    'datfile', {'a1.0000_b1.0000_d16.0.dat',   'a2.0000_b0.5000_d16.0.dat',   'a0.2500_b4.0000_d16.0.dat'} ...
);

baseDir  = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
phi_e    = linspace(0, 2*pi,    200);
phi_act  = linspace(pi/2, 3*pi/2, 200);

axs = gobjects(1, 3);
for k = 1:3
    axs(k) = axes(fig);
    hold(axs(k), 'on'); box(axs(k), 'on');
    axs(k).XTick = []; axs(k).YTick = [];
end

for k = 1:3
    c       = cases(k);
    ax      = axs(k);
    dataDir = fullfile(baseDir, c.dir);

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

    imagesc(ax, xplt, yplt, phi);
    ax.YDir = 'normal';
    axis(ax, 'equal'); axis(ax, 'tight');
    colormap(ax, cmocean('tempo'));
    clim(ax, [0 1]);

    B   = readmatrix(fullfile(dataDir, c.datfile));
    row = B(end-15, :);
    xcm = [row(2), row(5)];
    ycm = [row(3), row(6)];
    th  = [row(4), row(7)];

    ex0    = c.a * cos(phi_e);   ey0    = c.b * sin(phi_e);
    ex_act = c.a * cos(phi_act); ey_act = c.b * sin(phi_act);

    for p = 1:2
        ct = cos(th(p)); st = sin(th(p));
        patch(ax, ct*ex0 - st*ey0 + xcm(p), st*ex0 + ct*ey0 + ycm(p), ...
              'w', 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 2);
        plot(ax, ct*ex_act - st*ey_act + xcm(p), st*ex_act + ct*ey_act + ycm(p), ...
             '-', 'Color', col, 'LineWidth', 3);
    end

    before = ax.Children;
    streamslice(ax, xplt, yplt, U, V, 2);
    h = setdiff(ax.Children, before);
    set(h, 'Color', 'w', 'LineWidth', 1.5);
end

cb = colorbar(axs(3), 'Ticks', 0:0.5:1);
cb.Title.String = '\phi';

% --- Layout ------------------------------------------------------------------
ml = 0.3; mr = 0.9; mb = 0.3; mt = 0.3; gap = 0.2;
fig.Units = 'inches'; figW = fig.Position(3); figH = fig.Position(4);
set(axs, 'Units', 'inches');
panelH = figH - mb - mt;
panelW = (figW - ml - mr - 2*gap) / 3;

for k = 1:3
    axs(k).Position = [ml + (k-1)*(panelW + gap), mb, panelW, panelH];
end

cb.Units = 'inches';
cb.Position(1) = axs(3).Position(1) + axs(3).Position(3) + 0.1;

panel_label('\it (a)', axs(1), -0.1, 0.1);
panel_label('\it (b)', axs(2), -0.1, 0.1);
panel_label('\it (c)', axs(3), -0.1, 0.1);

function a = panel_label(str, ax, xOffset, yOffset)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', ax.FontSize, 'Units', 'inches');
    a.Position(1) = ax.Position(1) + xOffset;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + yOffset;
end
