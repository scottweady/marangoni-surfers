% fig_graphical_abstract.m — graphical abstract for the JFM table of contents: the unstable
% tandem chain of elongated surfers on the surfactant concentration, from fig_tandem, cropped
% to the journal's landscape aspect ratio with no axes or text. Written as a 300 dpi JPEG.
% The concentration is from tandem_snapshot.jl on a taller grid than fig_tandem's, so the
% window, which spans the full height of the grid, takes in more of the chain.

close all

ab_str = '4.0';
t_snap = 1000;
xc     = 37;          % center of the cropped window, in the frame of the chain centroid
sz     = [6 5];      % width and height, cm
outDir_ga = fullfile(fileparts(mfilename('fullpath')), 'output');

figsDir = fileparts(mfilename('fullpath'));
addpath(figsDir);
dataRoot = fullfile(figsDir, '..', 'data');
dataDir = fullfile(dataRoot, 'tandem');

% concentration field
fid = fopen(fullfile(dataRoot, 'graphical_abstract', sprintf('tandem_snapshot_%s_t%d.dat', ab_str, t_snap)));
hdr = fgetl(fid);
Nx  = str2double(regexp(hdr, 'Nx=(\d+)', 'tokens', 'once'));
Ny  = str2double(regexp(hdr, 'Ny=(\d+)', 'tokens', 'once'));
x   = str2num(fgetl(fid));  %#ok<ST2NM>
y   = str2num(fgetl(fid));  %#ok<ST2NM>
c   = fscanf(fid, '%f', [Nx, Ny])';
fclose(fid);

% particle configuration at the snapshot time; columns: t  x  y  θ  vx  vy  ω
dir_c  = fullfile(dataDir, ['tandem_' ab_str '_N50_d10.0_m5_ε0.01']);
f      = dir(fullfile(dir_c, 'frame_*.dat'));
frames = cellfun(@(f) load(fullfile(dir_c, f)), sort({f.name}), 'UniformOutput', false);
t      = cellfun(@(d) d(1,1), frames);
[~, n] = min(abs(t - t_snap));
d      = frames{n};

% same colormap as fig_tandem
c_rng = [0.4 0.6];
cmap  = cmocean('dense', 256);
cmap  = cmap(round(1 + 255 * min(max((linspace(0, 1, 1001) - c_rng(1)) / diff(c_rng), 0), 1)), :);

fig = figure('Units', 'centimeters', 'Position', [2 2 sz], 'Color', 'w', ...
             'PaperUnits', 'centimeters', 'PaperSize', sz, 'PaperPosition', [0 0 sz], 'InvertHardcopy', 'off');
ax  = axes(fig, 'Position', [0 0 1 1]);  hold(ax, 'on');
imagesc(ax, x, y, c, 'AlphaData', ~isnan(c));

a = sqrt(str2double(ab_str));  b = 1/a;
phi = linspace(0, 2*pi, 81);
for i = 1:size(d, 1)
    R  = [cos(d(i,4)) -sin(d(i,4)); sin(d(i,4)) cos(d(i,4))];
    xy = R * [a*cos(phi); b*sin(phi)];
    fill(ax, d(i,2) - mean(d(:,2)) + xy(1,:), d(i,3) - mean(d(:,3)) + xy(2,:), 'w', ...
         'EdgeColor', 'k', 'LineWidth', 0.4);
end

H = y(end) - y(1);  W = H * sz(1) / sz(2);
set(ax, 'YDir', 'normal', 'DataAspectRatio', [1 1 1], 'XLim', xc + W/2*[-1 1], 'YLim', y([1 end]));
axis(ax, 'off');
colormap(ax, cmap);
clim(ax, [0 1]);

jpgFile = fullfile(outDir_ga, 'fig_graphical_abstract.jpg');
print(fig, jpgFile, '-djpeg100', '-r300');

% print leaves the JFIF density unset; record the print resolution in dots per inch
fid = fopen(jpgFile, 'r+');
fseek(fid, 13, 'bof');
fwrite(fid, [1 1 44 1 44], 'uint8');
fclose(fid);
