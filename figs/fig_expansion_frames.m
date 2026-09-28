% fig_expansion_frames.m — render one PNG per saved timestep of the
% particle configuration shown in fig_expansion.m panel (a), for animation.

close all;

a = sqrt(1.0);
b = 1.0 / a;

% ── Load frames ───────────────────────────────────────────────────────────────

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'frames');
outDir  = fullfile(fileparts(mfilename('fullpath')), 'output', 'frames_expansion');
files   = dir(fullfile(dataDir, 'frame_*.dat'));
files   = sort({files.name});
Nframes = numel(files);

frames = cell(Nframes, 1);
for k = 1:Nframes
    frames{k} = load(fullfile(dataDir, files{k}));  % columns: t  x_re  x_im  θ  vx  vy  ω
end

mkdir(outDir);

% ── Fixed axes from final frame ───────────────────────────────────────────────

last  = frames{end};
pad   = 2 * max(a, b);
xlims = [min(last(:,2))-pad,  max(last(:,2))+pad];
ylims = [min(last(:,3))-pad,  max(last(:,3))+pad];

% ── Ellipse template (local frame) ───────────────────────────────────────────

s  = linspace(0, 2*pi, 64);
ex = a * cos(s);
ey = b * sin(s);

% ── Figure ───────────────────────────────────────────────────────────────────

if ~exist('P', 'var'), P = palette(); end
col = P.col;

fig = journal_figure([3.25 3.25], 2);
ax  = axes(fig); axis(ax, 'equal'); box(ax, 'on');
set(ax, 'XTick', [], 'YTick', []);

for k = 1:Nframes
    cla(ax); hold(ax, 'on');

    fdata = frames{k};
    x  = fdata(:,2);  y = fdata(:,3);  th = fdata(:,4);
    for i = 1:numel(x)
        R  = [cos(th(i)) -sin(th(i)); sin(th(i)) cos(th(i))];
        xy = R * [ex; ey];
        fill(ax, x(i)+xy(1,:), y(i)+xy(2,:), col, 'EdgeColor', 'none');
    end

    xlim(ax, xlims); ylim(ax, ylims);
    text(ax, 0.03, 0.97, sprintf('{\\it t} = %.2f', fdata(1,1)), ...
         'Units', 'normalized', 'HorizontalAlignment', 'left', ...
         'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'Margin', 2, ...
         'FontSize', ax.FontSize);

    exportgraphics(fig, fullfile(outDir, sprintf('frame_%04d.png', k)), 'Resolution', 300);
end

close(fig);
fprintf('Wrote %d frames to %s\n', Nframes, outDir);
