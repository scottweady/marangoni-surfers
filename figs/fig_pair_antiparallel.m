% fig_pair_antiparallel.m — anti-parallel ellipse pairs: (a)-(c) concentration + velocity
% fields at d0 = 8, (d) separation and (e) orientation vs time across aspect ratio.

close all;

P = palette();

fig = journal_figure([6.75, 4.35], 2, 'tex');

% ── Fields ────────────────────────────────────────────────────────────────────
dataDir = fullfile(dataRoot, 'antiparallel');
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
cb.Title.String = '{\it{c}}';


% ── Separation and orientation vs time ───────────────────────────────────────
traj_ab = [0.25, 1, 10];
traj_d  = [5.0, 4.0, 1.6];
Nt = numel(traj_ab);

% single-particle U, β, λ for the far-field ODEs (5.2)-(5.3)
E = readmatrix(fullfile(dataRoot, 'single', 'equilibrium_separation.dat'), 'NumHeaderLines', 1);

cmap = P.cmap(Nt+1);
cmap = cmap(2:end, :);
axSep = axes(fig); hold(axSep, 'on'); box(axSep, 'on');
axOmg = axes(fig); hold(axOmg, 'on'); box(axOmg, 'on');

for k = 1:Nt
    a_ell = sqrt(traj_ab(k));  b_ell = 1.0/a_ell;
    A    = readmatrix(fullfile(dataDir, sprintf('a%.4f_b%.4f_d%.1f.dat', a_ell, b_ell, traj_d(k))));
    A    = A(all(~isnan(A), 2), :);
    c    = cmap(k, :);
    t    = A(:,1);
    sep  = hypot(A(:,2) - A(:,5), A(:,3) - A(:,6));
    ang1 = unwrap(A(:,4));
    plot(axSep, t, sep, 'Color', c);
    plot(axOmg, t, ang1 + pi/2, 'Color', c, 'HandleVisibility', 'off');

    disp(['Processing trajectory ' num2str(k) ' of ' num2str(Nt)]);
    [~, i] = min(abs(E(:,1) - traj_ab(k)));
    uS = E(i,4);  beta = E(i,5);  lam = E(i,6);
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    [tf, Y] = ode45(@(~, y) farfield_rhs(y, uS, beta, lam), [0 50], A(1, 2:7)', opts);
    sep_ff = hypot(Y(:,1) - Y(:,4), Y(:,2) - Y(:,5));
    plot(axSep, tf, sep_ff, '--', 'Color', c);
    plot(axOmg, tf, unwrap(Y(:,3)) + pi/2, '--', 'Color', c, 'HandleVisibility', 'off');

end

linkaxes([axSep, axOmg], 'x');
xlabel(axSep, 'time, {\it t}');
xlabel(axOmg, 'time, {\it t}');
ylabel(axSep, '|{\it\bfx}_{cm,1} - {\it\bfx}_{cm,2}|');
ylabel(axOmg, '{\it\theta}_{cm,1} + \pi/2');

xlim(axSep, [0 50]);
xlim(axOmg, [0 50]);
ylim(axOmg, [-1.2 0.4]);

% ── Layout ────────────────────────────────────────────────────────────────────
topH = 4.2; botH = 4.8;

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
ml = 1.0; mr = 0.75; mb = 0.75; mt = 0.3; gapH = 1.2;
panelW = (figW - ml - mr - gapH) / 2;
panelH = botH - mb - mt;
axSep.Position = [ml,                 mb, panelW, panelH];
axOmg.Position = [ml + panelW + gapH, mb, panelW, panelH];

% ── Legend: simulation (solid) over far field (dashed), one label per a/b ────
axLeg = axes(fig, 'Units', 'inches', 'Position', axSep.Position, ...
             'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
hold(axLeg, 'on');
x0 = 0.04; L = 0.12; dy = 0.025; rowH = 0.1;
for k = 1:Nt
    yk = 0.93 - (k-1)*rowH;
    plot(axLeg, x0 + [0 L], yk + dy*[1 1], '-',  'Color', cmap(k,:));
    plot(axLeg, x0 + [0 L], yk - dy*[1 1], '--', 'Color', cmap(k,:));
    text(axLeg, x0 + L + 0.03, yk, ['{\it a/b} = ' num2str(traj_ab(k), '%.2f')], ...
         'VerticalAlignment', 'middle', 'FontName', axSep.FontName, 'FontSize', axSep.FontSize);
end

% ── Panel labels ──────────────────────────────────────────────────────────────
letters = 'abc';
for k = 1:Np
    a = aspan(k); b = 1.0/a;
    panel_label(sprintf('\\it (%s) {\\it a/b} \\rm = %.2f', letters(k), a/b), ax(k), -0.2, -0.1);
end
panel_label('\it (d)', axSep, -0.95, 0.2);
panel_label('\it (e)', axOmg, -0.95, 0.2);

% Far-field ODEs (5.2)-(5.3) for a pair, y = [x1 y1 θ1 x2 y2 θ2], r_ij = x_i - x_j.
function dy = farfield_rhs(y, uS, beta, lam)
    X = [y(1) y(2); y(4) y(5)];  th = [y(3) y(6)];
    dy = zeros(6, 1);
    for i = 1:2
        j  = 3 - i;
        p  = [cos(th(i)); sin(th(i))];  q = [-p(2); p(1)];
        pj = [cos(th(j)); sin(th(j))];
        r  = (X(i,:) - X(j,:))';  d = norm(r);  e = r/d;
        u  = uS*p + beta*(pj - e*(e'*pj))/(4*pi*d);
        w  = -beta*lam*( -3*(q'*e)*(e'*p)*(e'*pj) + (q'*e)*(pj'*p) + (q'*pj)*(e'*p) )/(4*pi*d^2);
        dy(3*i-2:3*i) = [u; w];
    end
end

function a = panel_label(str, ax, xOffset, yOffset)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', 20, 'Units', 'inches');
    a.Position(1) = ax.Position(1) + xOffset;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + yOffset;
end
