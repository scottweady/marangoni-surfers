% fig_schematic.m — 3D schematic of the Marangoni surfer setup:
% a fluid-filled box with the surfactant concentration field on the interface,
% ellipsoidal particles resting on top, and the bulk flow field below
% (examples/suspension/suspension_render.jl).

close all;

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'examples', 'suspension', 'output');
config  = 'lattice';

fid = fopen(fullfile(dataDir, sprintf('fields_%s.dat', config)));
hdr = fgetl(fid);
Nx  = str2double(regexp(hdr, 'Nx=(\d+)', 'tokens', 'once'));
Ny  = str2double(regexp(hdr, 'Ny=(\d+)', 'tokens', 'once'));
N   = str2double(regexp(hdr, 'N=(\d+)', 'tokens', 'once'));
a   = str2double(regexp(hdr, 'a=([\d.]+)', 'tokens', 'once'));
b   = str2double(regexp(hdr, 'b=([\d.]+)', 'tokens', 'once'));
alpha_arc = str2double(regexp(hdr, 'alpha=([\d.]+)', 'tokens', 'once'));

P = palette();

xplt = str2num(fgetl(fid)); %#ok<ST2NM>
yplt = str2num(fgetl(fid)); %#ok<ST2NM>
xcm  = str2num(fgetl(fid)); %#ok<ST2NM>
ycm  = str2num(fgetl(fid)); %#ok<ST2NM>
th   = str2num(fgetl(fid)); %#ok<ST2NM>

xL = 0.75 * min(xplt); xR = 0.75 * max(xplt);
yL = 0.75 * min(yplt); yR = 0.75 * max(yplt);
zB = 8; zT = 0;

phi = fscanf(fid, '%f', [Nx, Ny])';
fscanf(fid, '%f', [Nx, Ny]);   % u (unused)
fscanf(fid, '%f', [Nx, Ny]);   % v (unused)
fclose(fid);

Lz = zB;

fig = journal_figure([6.75 4], 2, 'tex');
ax  = axes(fig); hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');
ax.Position = [0 0 1 1];

draw_box(ax, xplt, yplt, Lz);

[X, Y] = meshgrid(xplt, yplt);
surf(ax, X, Y, zeros(size(phi)), phi, 'EdgeColor', 'none', 'FaceAlpha', 0.75);
colormap(ax, cmocean('tempo'));
clim(ax, [0 1]);

t   = linspace(0, 2 * pi, 256);
t_act = linspace(alpha_arc, 2 * pi - alpha_arc, 256);   % Dirichlet (active) arc
for k = 1:N
    ct = 1.1*cos(th(k)); st = 1.1*sin(th(k));

    ex = a * cos(t) * ct - b * sin(t) * st + xcm(k);
    ey = a * cos(t) * st + b * sin(t) * ct + ycm(k);
    ez = 0.01 * ones(size(t));   % lift above the surface — avoids z-fighting with surf()
    fill3(ax, ex, ey, ez, 'w', 'EdgeColor', 'k', 'LineWidth', 2);

    ex_act = a * cos(t_act) * ct - b * sin(t_act) * st + xcm(k);
    ey_act = a * cos(t_act) * st + b * sin(t_act) * ct + ycm(k);
    ez_act = 0.011 * ones(size(t_act));
    plot3(ax, ex_act, ey_act, ez_act, '-', 'Color', P.col, 'LineWidth', 3);
end

bfid = fopen(fullfile(dataDir, sprintf('bulk_%s.dat', config)));
bhdr = fgetl(bfid);
fclose(bfid);
NxB = str2double(regexp(bhdr, 'Nx=(\d+)', 'tokens', 'once'));
NyB = str2double(regexp(bhdr, 'Ny=(\d+)', 'tokens', 'once'));
NzB = str2double(regexp(bhdr, 'Nz=(\d+)', 'tokens', 'once'));

bulk = readmatrix(fullfile(dataDir, sprintf('bulk_%s.dat', config)), 'NumHeaderLines', 1);
grid3 = @(col) permute(reshape(bulk(:, col), [NxB, NyB, NzB]), [2 1 3]);
[Xb, Yb, Zb, Ub, Vb, Wb] = deal(grid3(1), grid3(2), grid3(3), grid3(4), grid3(5), grid3(6));

skip = 2;
ds = @(A) A(1:skip:end, 1:skip:end, 1:skip:end);
[Xb, Yb, Zb, Ub, Vb, Wb] = deal(ds(Xb), ds(Yb), ds(Zb), ds(Ub), ds(Vb), ds(Wb));

% keep = Zb >= -Lz;
% quiver3(ax, Xb(keep), Yb(keep), Zb(keep), Ub(keep), Vb(keep), Wb(keep), ...
    % 3, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);

% z-depth labels at the left corner of the box
x0 = min(xplt); y0 = min(yplt)-2;
tb1 = text(ax, x0, y0, 0+1.5, '{\it z} = 0', 'FontSize', 32, 'HorizontalAlignment', 'left');
tb2 = text(ax, x0, y0, -Lz-2, '{\it z} \rightarrow -\infty', 'FontSize', 32, 'HorizontalAlignment', 'left');
tb1.Position(1) = tb2.Position(1);

% Domain labels: S on the interface, V in the bulk
x1 = -2 + max(xplt); y1 = max(yplt);
text(ax, 8, 8, 0.5, '{\it S}', 'FontSize', 32, 'HorizontalAlignment', 'center');
text(ax, 14, 6, -0.2 * Lz, '{\it V}', 'FontSize', 32, 'HorizontalAlignment', 'center');

view(ax, [60 20]);
camzoom(ax, 1.6);
% xlim([xL xR])
% ylim([yL yR])


function draw_box(ax, xplt, yplt, Lz)
    x0 = min(xplt); x1 = max(xplt);
    y0 = min(yplt); y1 = max(yplt);

    cmap = cmocean('tempo');
    cmap = mean(cmap, 1);

    % Side walls (fluid), translucent — no solid bottom: the domain is
    % semi-infinite in depth, so the bottom is indicated with dotted lines only
    walls = {
        [x0 x1 x1 x0], [y0 y0 y0 y0], [0 0 -Lz -Lz]
        [x0 x1 x1 x0], [y1 y1 y1 y1], [0 0 -Lz -Lz]
        [x0 x0 x0 x0], [y0 y1 y1 y0], [0 0 -Lz -Lz]
        [x1 x1 x1 x1], [y0 y1 y1 y0], [0 0 -Lz -Lz]
        };
    for k = 1:size(walls, 1)
        patch(ax, walls{k, 1}, walls{k, 2}, walls{k, 3}, ...
            cmap, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    end

    % Dotted vertical edges and bottom rim — domain is semi-infinite in depth
    corners = [x0 y0; x1 y0; x1 y1; x0 y1];
    for k = 1:4
        plot3(ax, corners(k, 1) * [1 1], corners(k, 2) * [1 1], [0 -Lz], ':', ...
            'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    end

    rim = [corners; corners(1, :)];
    plot3(ax, rim(:, 1), rim(:, 2), -Lz * ones(5, 1), ':', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 2);

    plot3(ax, rim(:, 1), rim(:, 2), 0 * ones(5, 1), ':', ...
        'Color', [0.2 0.2 0.2], 'LineWidth', 2);
end
