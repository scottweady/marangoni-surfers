% fig_multipole_flows.m  —  Three-panel illustration of multipole expansion
%   (a) monopole flow: u = -r_hat
%   (b) dipole  flow: u = -(I - r_hat*r_hat)*e_y / r
%   (c) monopole coefficient vs separation (from sigma_vs_separation.dat)

close all

fig = journal_figure([6.75 2], 2, 'tex');

P = palette();

% ── Shared grid ───────────────────────────────────────────────────────────────

L  = 5;
Nq = 32;

[X, Y] = meshgrid(linspace(-L, L, Nq), linspace(-L, L, Nq));
R2 = X.^2 + Y.^2;  R = sqrt(R2);

rp   = 1.0;
mask = R < rp;

th   = linspace(0, 2*pi, 256);
cmap = P.cmap();
clim_max = 1/rp;

% ── Fine grid for pcolor ──────────────────────────────────────────────────────

Nc = 200;
[Xc, Yc] = meshgrid(linspace(-L, L, Nc), linspace(-L, L, Nc));
R2c = Xc.^2 + Yc.^2;  Rc = sqrt(R2c);
maskc = Rc < rp;

% ── (a) Monopole: u = -r_hat/r,  |u| = 1/r ───────────────────────────────────

Ua = -X ./ R;  Va = -Y ./ R;
Ua(mask) = NaN;  Va(mask) = NaN;

Mag_a = 0.125*ones(size(Rc));  Mag_a(maskc) = NaN;

sp1 = axes(fig);
hold(sp1, 'on'); box(sp1, 'on'); axis(sp1, 'equal');
colormap(sp1, cmap);  clim(sp1, [0, clim_max]);

contourf(sp1, Xc, Yc, Mag_a, 8, 'LineColor', 'none');
Ns_a    = 32;
r_arr_a = 2:6;   % radii at which to place arrows
for th_a = linspace(0, 2*pi*(1-1/Ns_a), Ns_a)
    r_vals = linspace(rp, L*sqrt(2), 300);
    xline = r_vals*cos(th_a);  yline = r_vals*sin(th_a);
    i_out = find(abs(xline) > L | abs(yline) > L, 1, 'first');
    if isempty(i_out); i_out = numel(r_vals)+1; end
    plot(sp1, xline(1:i_out-1), yline(1:i_out-1), 'Color', [0.15 0.15 0.15], 'LineWidth', 1);
    for r_a = r_arr_a
        if r_a < r_vals(i_out-1)
            arrowhead(sp1, r_a*cos(th_a), r_a*sin(th_a), pi+th_a);
        end
    end
end

fill(sp1, rp*cos(th), rp*sin(th), [0.88 0.88 0.88], 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 2);
phi_act = linspace(pi/2, 3*pi/2, 200);
plot(sp1, -rp*sin(phi_act), rp*cos(phi_act), '-', 'Color', P.col, 'LineWidth', 3);

axis(sp1, [-L L -L L]);
xlabel(sp1, '{\it x}');
ylabel(sp1, '{\it y}');

% ── (b) Dipole: u = -(I - r_hat*r_hat)*e_y / r ───────────────────────────────
%   u_x = x*y/r^3,  u_y = -x^2/r^3,  |u| = |x|/r^2

R3 = R2 .* R;
Ub = X .* Y ./ R3;  Vb = -X.^2 ./ R3;
Ub(mask) = NaN;  Vb(mask) = NaN;

Mag_b = abs(Xc) ./ R2c;  Mag_b(maskc) = NaN;

sp2 = axes(fig);
hold(sp2, 'on'); box(sp2, 'on'); axis(sp2, 'equal');
colormap(sp2, cmap);  clim(sp2, [0, clim_max]);

% pcolor(sp2, Xc, Yc, Mag_b); shading(sp2, 'interp');
contourf(sp2, Xc, Yc, Mag_b, 8, 'LineColor', 'none');
th_c    = linspace(0, 2*pi, 1000);
dth = pi/16 * (1:2:17);
th_arr_b = [pi/2 + dth, pi/2 - dth];  % angles at which to place arrows

for r_c = linspace(2, sqrt(2) * L, 8)
    xc = r_c*cos(th_c);  yc = r_c*sin(th_c);
    in_box = abs(xc) <= L & abs(yc) <= L;
    xc(~in_box) = NaN;  yc(~in_box) = NaN;
    plot(sp2, xc, yc, 'Color', [0.15 0.15 0.15], 'LineWidth', 1);
    for th_b = th_arr_b
        xa = r_c*cos(th_b);  ya = r_c*sin(th_b);
        if abs(xa) > L || abs(ya) > L; continue; end
        arrowhead(sp2, xa, ya, atan2(-(1+cos(2*th_b)), sin(2*th_b)));
    end
end

fill(sp2, rp*cos(th), rp*sin(th), [0.88 0.88 0.88], 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 2);
plot(sp2, -rp*sin(phi_act), rp*cos(phi_act), '-', 'Color', P.col, 'LineWidth', 3);

axis(sp2, [-L L -L L]);
xlabel(sp2, '{\it x}');
set(sp2, 'YTickLabel', {});

cb = colorbar(sp2, 'Ticks', [0, clim_max], 'TickLabels', {'0', '|{\it\bf u}|_{max}'});

% ── (c) Monopole coefficient vs separation ────────────────────────────────────

dataDir = fullfile(dataRoot, 'multipole');
data     = load(fullfile(dataDir, 'sigma_vs_separation.dat'));
d        = data(:,1);
I1_th0   = data(:,2);
I1_thpi2 = data(:,4);
d_fine   = linspace(d(1), d(end), 300)';

fit0 = 1./(d_fine.*log(d_fine));
fit0 = fit0 * abs(I1_th0(end)/fit0(end));
sp3 = axes(fig);
hold(sp3, 'on'); box(sp3, 'on');

plot(sp3, d, abs(I1_th0),   'Color', P.c1, 'DisplayName', '{\it \theta} = 0');
plot(sp3, d, abs(I1_thpi2), 'Color', P.c2, 'DisplayName', '{\it \theta} = {\it \pi}/2');
plot(sp3, d_fine, fit0, '--', 'Color', [0.6 0.2 0.2], 'DisplayName', '\propto 1/({\it{r}} log{\it r})')
xlabel(sp3, 'separation, {\it r}');
ylabel(sp3, ['{\it{m}}_{0,1}']);
legend(sp3, 'Location', 'northeast', 'Color', 'none', 'FontSize', sp3.FontSize);
xlim(sp3, [0 100]);
ylim(sp3, [-0.02 0.2]);

% ── Layout ────────────────────────────────────────────────────────────────────

for ax = [sp1, sp2, sp3]
    ax.Units = 'inches';
    ax.ActivePositionProperty = 'position';
end
drawnow;

lm  = 0.55;   % left margin
rm  = 0.25;   % right margin
bm  = 0.85;   % bottom margin
tm  = 0.25;   % top margin
cgL = 0.4;   % gap between flow panels
cgR = 1.5;   % gap before data panel

figW = fig.PaperSize(1);
figH = fig.PaperSize(2);

panelH  = figH - bm - tm;
flowW   = panelH;                            % square flow panels
dataW   = figW - lm - rm - 2*flowW - cgL - cgR;

sp1.Position = [lm,                          bm, flowW, panelH];
sp2.Position = [lm + flowW + cgL,            bm, flowW, panelH];
sp3.Position = [lm + 2*flowW + cgL + cgR,   bm, dataW,  panelH];

drawnow;
sp1.Position = [lm,                          bm, flowW, panelH];
sp2.Position = [lm + flowW + cgL,            bm, flowW, panelH];
sp3.Position = [lm + 2*flowW + cgL + cgR,   bm, dataW,  panelH];

cb.Units = 'inches';
cb.Position = [sp2.Position(1) + sp2.Position(3) + 0.12, bm, 0.12, panelH];

% ── Panel labels ──────────────────────────────────────────────────────────────

a = panel_label('{\it (a)}', sp1, -0.65, 0.05);
b = panel_label('{\it (b)}', sp2, -0.5, 0.05);
c = panel_label('{\it (c)}', sp3, -0.11, 0.05);

function arrowhead(ax, x, y, ang)
    sz = 0.25;
    hw = sz * 0.5;
    vx = sz*[-1, -1, 0];  vy = hw*[-1, 1, 0];
    patch(ax, vx*cos(ang) - vy*sin(ang) + x, ...
              vx*sin(ang) + vy*cos(ang) + y, ...
              [0.15 0.15 0.15], 'EdgeColor', 'none');
end

function a = panel_label(str, ax, xOffset, yOffset)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', ax.FontSize, 'Units', 'inches');
    a.Position(1) = ax.Position(1) + xOffset;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + yOffset;
end
