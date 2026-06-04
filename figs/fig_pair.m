% Combined figure:
%   Panel (a): surfactant concentration field for a=2 pair (midpoint snapshot)
%   Panel (b): boundary-to-boundary distance for a=2,3,4
%
% Requires: concentration_snapshot.dat (from concentration_snapshot.jl)
%           fig3/fig3_a*_b*_d*.dat

clear; clc; close all;

baseDir = fileparts(mfilename('fullpath'));

% ---- Shared colormap: 4 colors for a = 1,2,3,4 -----------------------------
cmap = flipud(cmocean('tempo', 5));
cmap = cmap(1:end-1,:);

% ---- Figure -----------------------------------------------------------------
fig = journal_figure([6.5, 2.25], 2, 'tex');

axConc = subplot(1, 2, 1); hold(axConc, 'on'); box(axConc, 'on');
axConc.XColor = 'k'; axConc.YColor = 'k';

axSep  = subplot(1, 2, 2); hold(axSep,  'on'); box(axSep,  'on');

% ==========================================================================
% Panel (a): concentration field
% ==========================================================================
snapFile = fullfile(baseDir, 'concentration_snapshot.dat');
fid = fopen(snapFile);
if fid == -1
    error('Cannot open %s — run concentration_snapshot.jl first.', snapFile);
end
hdr  = fgetl(fid);
Nx   = str2double(regexp(hdr, 'Nx=(\d+)', 'tokens', 'once'));
Ny   = str2double(regexp(hdr, 'Ny=(\d+)', 'tokens', 'once'));
xplt = str2num(fgetl(fid));   %#ok<ST2NM>
yplt = str2num(fgetl(fid));   %#ok<ST2NM>
phi  = fscanf(fid, '%f', [Nx, Ny])';   % Ny x Nx
U    = fscanf(fid, '%f', [Nx, Ny])';
V    = fscanf(fid, '%f', [Nx, Ny])';
fclose(fid);

imagesc(axConc, xplt, yplt, phi);
axConc.YDir = 'normal';
axis(axConc, 'equal'); axis(axConc, 'tight');
colormap(axConc, cmocean('tempo'));
axConc.XTick = []; axConc.YTick = [];

% Overlay ellipse outlines at the snapshot configuration
a_ell = 2.0; b_ell = 0.5;
trajFile = fullfile(baseDir, 'fig3', 'fig3_a2.0_b0.5000_d16.0.dat');
B   = readmatrix(trajFile);
row = B(floor(size(B,1)/2), :);
xcm_snap = [row(2), row(5)];
ycm_snap = [row(3), row(6)];
th_snap  = [row(4), row(7)];

phi_e = linspace(0, 2*pi, 200);
ex0   = a_ell*cos(phi_e);
ey0   = b_ell*sin(phi_e);
for k = 1:2
    ct = cos(th_snap(k)); st = sin(th_snap(k));
    xe = ct*ex0 - st*ey0 + xcm_snap(k);
    ye = st*ex0 + ct*ey0 + ycm_snap(k);
    patch(axConc, xe, ye, 'w', 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 2);
end

% Surface velocity streamlines
h = streamslice(axConc, xplt, yplt, U, V, 1);
set(h, 'Color', 'w', 'LineWidth', 1.5);

cb = colorbar(axConc, 'Location', 'westoutside', 'Ticks', 0:0.5:1)
cb_label = ylabel(cb, '{\it \phi}', 'FontName', 'Times', 'FontSize', axConc.FontSize, 'Rotation', 0);
clim(axConc, [0 1]);

% ==========================================================================
% Panel (b): separation vs time for a=2,3,4
% ==========================================================================
dataDir_ellipse = fullfile(baseDir, 'fig3');
files_e = dir(fullfile(dataDir_ellipse, 'fig3_a*_b*_d*.dat'));
aVals = nan(numel(files_e),1); bVals = nan(numel(files_e),1); match = false(numel(files_e),1);
for k = 1:numel(files_e)
    dVal = parse_param(files_e(k).name, 'd');
    aVals(k) = parse_param(files_e(k).name, 'a');
    bVals(k) = parse_param(files_e(k).name, 'b');
    if ~isnan(dVal); match(k) = true; end
end

filesD = files_e(match);
aValsD = aVals(match); bValsD = bVals(match);
[aValsD, idx] = sort(aValsD);
filesD = filesD(idx); bValsD = bValsD(idx);

keep   = arrayfun(@(a) any(abs(a - [2,3,4]) < 1e-9), aValsD);
filesD = filesD(keep); aValsD = aValsD(keep); bValsD = bValsD(keep);

if isempty(filesD)
    warning('Panel (b): no files found for a=2,3,4.');
else
    nGroups       = numel(filesD);
    sepHandles    = gobjects(nGroups, 1);
    legendEntries = strings(nGroups, 1);

    for k = 1:nGroups
        A = readmatrix(fullfile(filesD(k).folder, filesD(k).name));
        if isempty(A) || size(A,2) < 8; continue; end
        a = aValsD(k); b = bValsD(k);
        sepHandles(k)    = plot(axSep, A(:,1), A(:,8), 'Color', cmap(round(a),:));
        legendEntries(k) = strcat('{\it a/b} = ', num2str(round(a/b)));
    end

    xlabel(axSep, 'time, {\it t}');
    ylabel(axSep, 'dist({\it \Gamma}_1, {\it \Gamma}_2)');
    valid = isgraphics(sepHandles) & strlength(legendEntries) > 0;
    legend(axSep, sepHandles(valid), cellstr(legendEntries(valid)), ...
           'Location', 'northeast', 'Color', 'none');
end

% ==========================================================================
% Layout
% ==========================================================================
axConc.ActivePositionProperty = 'position';
axSep.ActivePositionProperty  = 'position';
axConc.Units = 'inches'; axSep.Units = 'inches';

tiS = axSep.TightInset;

leftMargin   = max(0.8,  tiS(1)) + 0.2;
rightMargin  = max(0.22, tiS(3)) + 0.06;
bottomMargin = max(0.42, tiS(2)) + 0.08;
topMargin    = max(0.24, tiS(4)) + 0.06;
midGap       = 1.2;

panelW = (fig.PaperSize(1) - leftMargin - rightMargin - midGap) / 2;
panelH = fig.PaperSize(2) - bottomMargin - topMargin;

axConc.Position = [leftMargin,                   bottomMargin, panelW, panelH];
axSep.Position  = [leftMargin + panelW + midGap, bottomMargin, panelW, panelH];

axConc.Position = [leftMargin,                   bottomMargin, panelW, panelH];
axSep.Position  = [leftMargin + panelW + midGap, bottomMargin, panelW, panelH];

cb.Units = 'inches';
cb_label.Units = 'inches';
cb_label.Position(1) = cb.Position(1) - 1.4;
% cb_label.Position(2) = 0.5*(cb.Position(2) + cb.Position(4));
cb_label.FontSize = axConc.FontSize;

% Panel labels.
panel_label('\it (a)', axConc);
panel_label('\it (b)', axSep);

% ---- helpers ---------------------------------------------------------------
function val = parse_param(fname, key)
    token = regexp(fname, ['_' key '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'], 'tokens', 'once');
    if isempty(token); val = NaN; else; val = str2double(token{1}); end
end

function a = panel_label(str, ax)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', 18, 'Units', 'inches');
    a.Position(1) = ax.Position(1) - 0.9;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + 0.2;
end
