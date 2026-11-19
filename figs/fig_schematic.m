% fig_schematic.m — schematic of the model: (a) a suspension of surfers on the surface S
% of a semi-infinite fluid V, colored by the surfactant concentration from a simulation
% snapshot (schematic_snapshot.jl); (b) zoom on one particle showing the concentration on
% its boundary, its center of mass, and its orientation.

close all

% ── Load snapshot ────────────────────────────────────────────────────────────

dataDir = fullfile(dataRoot, 'schematic');
fid = fopen(fullfile(dataDir, 'schematic_snapshot.dat'));
hdr = fgetl(fid);
hv  = @(name) str2double(regexp(hdr, [' ' name '=([-\d.]+)'], 'tokens', 'once'));
Nx = hv('Nx'); Ny = hv('Ny'); N = hv('N'); Nz = hv('Nz'); Nb = hv('Nb');
a  = hv('a');  b  = hv('b');  k = hv('k');

xplt = fscanf(fid, '%f', [1 Nx]);
yplt = fscanf(fid, '%f', [1 Ny]);
xcm  = fscanf(fid, '%f', [1 N]);
ycm  = fscanf(fid, '%f', [1 N]);
th   = fscanf(fid, '%f', [1 N]);
phi  = fscanf(fid, '%f', [Nx Ny])';
xz   = fscanf(fid, '%f', [1 Nz]);
yz   = fscanf(fid, '%f', [1 Nz]);
phiz = fscanf(fid, '%f', [Nz Nz])';
sb   = fscanf(fid, '%f', [1 Nb]);
phib = fscanf(fid, '%f', [1 Nb]);
fclose(fid);

P = palette();
ellipse_xy = @(s, j) deal(xcm(j) + a*cos(s)*cos(th(j)) - b*sin(s)*sin(th(j)), ...
                          ycm(j) + a*cos(s)*sin(th(j)) + b*sin(s)*cos(th(j)));
s = linspace(0, 2*pi, 128);

fig = journal_figure([6.5 2.6], 2, 'tex');
set(fig, 'defaultTextFontSize', 20);

% ── (a) Suspension on the surface of a semi-infinite fluid ───────────────────

axA = axes(fig); hold(axA, 'on'); axis(axA, 'equal', 'off');

x0 = xplt(1); x1 = xplt(end); y0 = yplt(1); y1 = yplt(end);
Lz = 0.5 * (x1 - x0);
wall = mean(P.cmap(), 1);
walls = {[x0 x1 x1 x0], [y0 y0 y0 y0], [0 0 -Lz -Lz]
         [x0 x1 x1 x0], [y1 y1 y1 y1], [0 0 -Lz -Lz]
         [x0 x0 x0 x0], [y0 y1 y1 y0], [0 0 -Lz -Lz]
         [x1 x1 x1 x1], [y0 y1 y1 y0], [0 0 -Lz -Lz]};
for j = 1:size(walls, 1)
    patch(axA, walls{j, :}, wall, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
corners = [x0 y0; x1 y0; x1 y1; x0 y1; x0 y0];
for j = 1:4
    plot3(axA, corners(j, 1)*[1 1], corners(j, 2)*[1 1], [0 -Lz], ':', 'Color', [0.5 0.5 0.5]);
end
plot3(axA, corners(:, 1), corners(:, 2), -Lz*ones(5, 1), ':', 'Color', [0.5 0.5 0.5]);

surf(axA, [x0 x1; x0 x1], [y0 y0; y1 y1], zeros(2), phi, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
plot3(axA, corners(:, 1), corners(:, 2), zeros(5, 1), '-', 'Color', P.ink);

% in dark mode the particles are dark with a light edge, as in the presentation's animations
face = 'w'; edge = 'k';
if P.dark, face = [42 34 48] / 255; edge = P.fg; end
for j = 1:N
    [ex, ey] = ellipse_xy(s, j);
    fill3(axA, ex, ey, 0*ex, face, 'EdgeColor', edge, 'LineWidth', 1);
end

% Zoom window for (b)
hz = (xz(end) - xz(1)) / 2;
zx = xcm(k) + hz*[-1 1 1 -1 -1];
zy = ycm(k) + hz*[-1 -1 1 1 -1];
plot3(axA, zx, zy, zeros(1, 5), 'k-', 'LineWidth', 1.5);

text(axA, x1+2, y1-2, -3, '  {\it S}', 'VerticalAlignment', 'bottom');
text(axA, x1-2, y0, -0.85*Lz, '{\it V}    ', 'HorizontalAlignment', 'right');
text(axA, x0, y1, 0, '{\it z} = 0  ', 'HorizontalAlignment', 'right');
text(axA, x0, y1, -Lz, '{\it z} \rightarrow -\infty  ', 'HorizontalAlignment', 'right');

colormap(axA, P.cmap()); clim(axA, [0 1]);
view(axA, [-30 35]);

% ── (b) Zoom on particle k ───────────────────────────────────────────────────

axB = axes(fig); hold(axB, 'on'); axis(axB, 'equal'); box(axB, 'on');
set(axB, 'XTick', [], 'YTick', [], 'LineWidth', 1.5, 'Layer', 'top');

imagesc(axB, xz, yz, phiz, 'AlphaData', ~isnan(phiz));
axB.YDir = 'normal';
for j = 1:N
    [ex, ey] = ellipse_xy(s, j);
    fill(axB, ex, ey, 'w', 'EdgeColor', 'k', 'LineWidth', 1);
end

R     = [cos(th(k)) -sin(th(k)); sin(th(k)) cos(th(k))];
body  = @(s) R * [a*cos(s); b*sin(s)] + [xcm(k); ycm(k)];
tang  = @(s) R * [-a*sin(s); b*cos(s)];
normal = @(s) R * [b*cos(s); a*sin(s)];

% Dirichlet arc: thick black arc
pD = body(linspace(pi/2, 3*pi/2, 256));
plot(axB, pD(1,:), pD(2,:), 'k-', 'LineWidth', 2.5);

% Neumann arc: impermeable wall, hatched with inward normal strokes of length w at equal arc length
w  = 0.05;
sf = linspace(-pi/2, pi/2, 2000);
ds = vecnorm(tang(sf)) * (sf(2) - sf(1));
arclen = [0 cumsum(ds(1:end-1))];
sh = interp1(arclen, sf, 0.06:0.06:arclen(end));
nh = normal(sh) ./ vecnorm(normal(sh));
ph = body(sh);
qh = ph - w*nh;
plot(axB, reshape([ph(1,:); qh(1,:); nan(size(sh))], 1, []), ...
          reshape([ph(2,:); qh(2,:); nan(size(sh))], 1, []), 'k-', 'LineWidth', 0.75);

[ex, ey] = ellipse_xy(pi, k);
text(axB, xcm(k) - 0.1*(ex - xcm(k)), ycm(k) + 1.3*(ey - ycm(k)), '{\it\Gamma}_{\it{i,D}}', ...
     'HorizontalAlignment', 'center');
[ex, ey] = ellipse_xy(-pi/3, k);
text(axB, xcm(k) + 1.6*(ex - xcm(k)), ycm(k) + 1.4*(ey - ycm(k)), '{\it\Gamma}_{\it{i,N}}', ...
     'HorizontalAlignment', 'center');
text(axB, xcm(k) - 0.42, ycm(k) - 0.05, '\Omega_{\it{i}}', 'HorizontalAlignment', 'center');

% Center of mass, orientation p_i, and angle θ_cm,i
p  = [cos(th(k)) sin(th(k))];
Lp = 1.4 * a;
Lt = 1.6 * b;
plot(axB, xcm(k) + [0 Lp], ycm(k) + [0 0], 'k:', 'LineWidth', 2);
ta = linspace(0, th(k), 64);
plot(axB, xcm(k) + 0.3*a*cos(ta), ycm(k) + 0.3*a*sin(ta), 'k-', 'LineWidth', 1);
text(axB, xcm(k) + 0.5*a*cos(th(k)/2), ycm(k) + 0.5*a*sin(th(k)/2), '\theta_{cm,{\it i}}', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
quiver(axB, xcm(k), ycm(k), Lp*p(1), Lp*p(2), 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.2);
text(axB, xcm(k) + Lp*p(1), ycm(k) + Lp*p(2), '  {\bf p}_{\it i}');
plot(axB, xcm(k), ycm(k), 'k.', 'MarkerSize', 24);
text(axB, xcm(k) - 0.05, ycm(k) - 0.1, '{\bf x}_{cm,{\it i}}', 'VerticalAlignment', 'top');

% labels sit on the concentration field, so stay dark in either mode
set(findobj(axB, 'Type', 'text'), 'Color', 'k');

xlim(axB, xz([1 end])); ylim(axB, yz([1 end]));
colormap(axB, P.cmap()); clim(axB, [0 1]);
cb = colorbar(axB, 'Ticks', 0:0.5:1);

% ── Layout ───────────────────────────────────────────────────────────────────

figW = fig.PaperSize(1); figH = fig.PaperSize(2);
for h = [axA, axB, cb]
    h.Units = 'inches';
end
sz = 0.86*figH;
axB.Position = [figW - sz - 1.1, 0.07*figH, sz, sz];
cb.Position  = [figW - 0.95, 0.07*figH, 0.18, sz];
cb.Title.String = '{\it c}';
cb.Title.FontSize = axB.FontSize;
axA.Position = [0.1, -0.45, axB.Position(1) - 0.3, figH + 1.0];
camzoom(axA, 1.6);

annotation(fig, 'textbox', [0 0.9 0.1 0.1], 'String', '\it(a)', 'EdgeColor', 'none', 'FontSize', 20, 'FontName', 'Times');
annotation(fig, 'textbox', [(axB.Position(1) - 0.45)/figW 0.9 0.1 0.1], 'String', '\it(b)', 'EdgeColor', 'none', 'FontSize', 20, 'FontName', 'Times');
