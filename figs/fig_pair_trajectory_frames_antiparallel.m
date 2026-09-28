% fig_pair_trajectory_frames_antiparallel.m — render one PNG per saved
% timestep of a two-ellipse antiparallel trajectory (examples/pair/pair.jl
% output), for the aspect ratios shown in fig_pair_fields_antiparallel.

close all;

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'oblate');
outDir  = fullfile(fileparts(mfilename('fullpath')), 'output', 'frames_oblate');
aspan     = [0.25];
Delta     = 12.0;

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'parallel');
outDir  = fullfile(fileparts(mfilename('fullpath')), 'output', 'frames_parallel');
aspan     = [2.0];
Delta     = 16.0;

alpha_arc = pi/2;

P   = palette();
col = P.col;

phi_e   = linspace(0, 2*pi, 200);
phi_act = linspace(alpha_arc, 2*pi - alpha_arc, 200);

for a = aspan
    b = 1.0 / a;
    fname = sprintf('a%.4f_b%.4f_d%.1f.dat', a, b, Delta);
    A = readmatrix(fullfile(dataDir, fname));
    A = A(all(~isnan(A), 2), :);

    t   = A(:,1);
    x1  = A(:,2); y1 = A(:,3); th1 = A(:,4);
    x2  = A(:,5); y2 = A(:,6); th2 = A(:,7);

    ex0    = a*cos(phi_e);   ey0    = b*sin(phi_e);
    ex_act = a*cos(phi_act); ey_act = b*sin(phi_act);

    lim = 2 * (Delta/2 + a);

    caseDir = fullfile(outDir, sprintf('a%.4f_b%.4f_d%.1f', a, b, Delta));
    mkdir(caseDir);

    % fig = figure('Visible', 'off', 'Color', 'w');
    fig = journal_figure([2 1], 2);
    ax  = axes(fig); axis(ax, 'equal'); box(ax, 'on');

    for k = 1:1:min(numel(t), 200)
        cla(ax); hold(ax, 'on');
        ycm = [-x1(k), -x2(k)]; xcm = [y1(k), y2(k)]; th = [th1(k)-pi/2, th2(k)-pi/2];
        % xcm = [x1(k), x2(k)]; ycm = [y1(k), y2(k)]; th = [th1(k), th2(k)];
        for j = 1:2
            ct = cos(th(j)); st = sin(th(j));
            patch(ax, ct*ex0 - st*ey0 + xcm(j), st*ex0 + ct*ey0 + ycm(j), ...
                  'w', 'FaceColor', 0.9*[1 1 1], 'EdgeColor', 'k', 'LineWidth', 2);
            plot(ax, ct*ex_act - st*ey_act + xcm(j), st*ex_act + ct*ey_act + ycm(j), ...
                 '-', 'Color', col, 'LineWidth', 2);
        end
        ax.YLim = [-10, 10]; ax.XLim = [-5 55];
        % ax.YLim = [-15, 15]; ax.XLim = [-15 15];
        % text(ax, -lim*0.95, lim*0.4, sprintf('{\\it t} = %.1f', t(k)), 'FontName', 'Times');

        exportgraphics(fig, fullfile(caseDir, sprintf('frame_%04d.png', k)), 'Resolution', 300);
        
    end

    close(fig);
    fprintf('Wrote %d frames to %s\n', numel(t), caseDir);
end
