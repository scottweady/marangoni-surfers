% Combined figure:
%   Panel (a): circular-case (a=1) separation with sqrt fit (from plot_fig3_circle_fit)
%   Panel (b): boundary-to-boundary distance for a=2,3,4 (from plot_fig3_trajectory panel a)
%
% Expected filename pattern: fig3_a{a}_b{b}_d{d}.dat
% Expected columns: 1:t, 2:x1, 3:y1, 4:theta1, 5:x2, 6:y2, 7:theta2, 8:separation

clear; clc; close all;

baseDir = fileparts(mfilename('fullpath'));

% ---- Panel (a) settings: circular case ------------------------------------
dataDirs = {
    fullfile(baseDir, 'fig3')
};
dataLabels = {'results'};

% ---- Panel (b) settings: aspect-ratio comparison --------------------------
targetD_ellipse = 10.0;
dataDir_ellipse = fullfile(baseDir, 'fig3');

% ---- Shared colormap: 4 colors for a = 1,2,3,4 ---------------------------
cmap = (cmocean('ice', 5)).^0.5;
cmap = cmap(1:end-1,:);

% ---- Figure ---------------------------------------------------------------
fig = journal_figure([6.5, 2.25], 2, 'tex');

axCirc  = subplot(1, 2, 1); hold(axCirc,  'on'); box(axCirc,  'on');
axSep   = subplot(1, 2, 2); hold(axSep,   'on'); box(axSep,   'on');

axInsetCirc = axes('Parent', fig);
hold(axInsetCirc, 'on'); box(axInsetCirc, 'on');
axInsetCirc.Color = 'none';
axInsetCirc.LooseInset = [0, 0, 0, 0];

axInsetSep = axes('Parent', fig);
hold(axInsetSep, 'on'); box(axInsetSep, 'off');
axInsetSep.Color = 'none'; axInsetSep.XColor = 'none'; axInsetSep.YColor = 'none';
axInsetSep.LooseInset = [0, 0, 0, 0];

% ==========================================================================
% Panel (a): circular case
% ==========================================================================
filesUse = struct('folder', {}, 'name', {}, 'sourceLabel', {});
for idir = 1:numel(dataDirs)
    files = dir(fullfile(dataDirs{idir}, 'fig3_a*_b*_d*.dat'));
    for k = 1:numel(files)
        aVal = parse_param(files(k).name, 'a');
        dVal = parse_param(files(k).name, 'd');
        if ~isnan(aVal) && ~isnan(dVal) && abs(aVal - 1.0) < 1e-9
            filesUse(end+1) = struct( ...
                'folder', files(k).folder, ...
                'name',   files(k).name, ...
                'sourceLabel', dataLabels{idir}); %#ok<SAGROW>
        end
    end
end

  for k = 1:numel(filesUse)
      filePath = fullfile(filesUse(k).folder, filesUse(k).name);
      A = readmatrix(filePath);
      if isempty(A) || size(A,2) < 7; continue; end

      t  = A(:,1);
      x1 = A(:,2); y1 = A(:,3); th1 = A(:,4);
      x2 = A(:,5); y2 = A(:,6); th2 = A(:,7);
      sep = A(:,8);

      if numel(t) < 2; continue; end

      cLine = cmap(1,:);  % a=1
      plot(axCirc, t, sep, '-', 'Color', cLine, 'DisplayName', '{\it a/b} = 1');
      tmp = 1 - (sep / sep(1)).^2;
      id = 1 : floor(numel(tmp)/4);
      p = polyfit(t(id), tmp(id), 1);
      tc = 1/p(1);
      t_sqrt = linspace(0, tc, 100);
      sep_sqrt = 14.0 * sqrt(1 - t_sqrt/tc);
      plot(axCirc, t_sqrt, sep_sqrt, '--', 'Color', cLine, 'DisplayName', '\propto (1 - {\it t/t_{sqrt}})^{1/2}');

      id = (length(t) - 40) : length(t);
      p = polyfit(t(id), sep(id), 1);
      tf = -p(2)/p(1);
      [~,i0] = min(abs(t(id)-(tf-100)));
      t_lin = exp(linspace(log(tf-5e2), log(tf), 100));
      sep_lin = p(1)*t_lin + p(2);
      plot(axInsetCirc, tf-t, sep, '-', 'Color', cLine)
      plot(axInsetCirc, tf-t_lin, sep_lin, ':', 'Color', [0.2 0.2 0.2])
      xlabel(axInsetCirc, '{\it t_c} - {\it t}');

      % plot(axCirc, t_lin, sep_lin, ':', 'Color', [0.2 0.2 0.2], 'DisplayName', '\propto {\it t_{lin}} - {\it t}');
      
      axInsetCirc.XLim = [1e-1, 100];
      axInsetCirc.YLim = [1e-3, 1];
      axInsetCirc.XScale = 'log';
      axInsetCirc.YScale = 'log';

      % set(axInsetCirc, 'XScale', 'log', 'YScale', 'log');
      % aVal = 1.0; bVal = 1.0;
      % nSnap  = min(10, numel(t));
      % tSnap  = linspace(t(1), t(end), nSnap);
      % x1s = interp1(t, x1, tSnap); y1s = interp1(t, y1, tSnap);
      % x2s = interp1(t, x2, tSnap); y2s = interp1(t, y2, tSnap);
      % th1s = interp1(t, unwrap(th1), tSnap);
      % th2s = interp1(t, unwrap(th2), tSnap);
      % phi  = linspace(0, 2*pi, 80);
      % ex0  = cos(phi); ey0 = sin(phi);

      % for j = 1:nSnap
      %     alpha = 0.15 + 0.85*(j-1)/max(1, nSnap-1);
      %     cFill = blend_with_white(blend_with_white(cLine, alpha), 0.9);
      %     ct1 = cos(th1s(j)); st1 = sin(th1s(j));
      %     patch(axInsetCirc, x1s(j)+ct1*ex0-st1*ey0, y1s(j)+st1*ex0+ct1*ey0, ...
      %         cFill, 'FaceColor', cFill, 'EdgeColor', 'k', 'LineWidth', 0.7);
      %     ct2 = cos(th2s(j)); st2 = sin(th2s(j));
      %     patch(axInsetCirc, x2s(j)+ct2*ex0-st2*ey0, y2s(j)+st2*ex0+ct2*ey0, ...
      %         cFill, 'FaceColor', cFill, 'EdgeColor', 'k', 'LineWidth', 0.7);
      % end
  end

xlabel(axCirc, 'time, {\it t}');
ylabel(axCirc, 'dist({\it \Gamma}_1, {\it \Gamma}_2)');
legend(axCirc, 'Location', 'northeast', 'Color', 'none');

% ==========================================================================
% Panel (b): ellipse aspect-ratio comparison
% ==========================================================================
files_e = dir(fullfile(dataDir_ellipse, 'fig3_a*_b*_d*.dat'));
aVals = nan(numel(files_e),1); bVals = nan(numel(files_e),1); match = false(numel(files_e),1);
for k = 1:numel(files_e)
    dVal = parse_param(files_e(k).name, 'd');
    aVals(k) = parse_param(files_e(k).name, 'a');
    bVals(k) = parse_param(files_e(k).name, 'b');
    if ~isnan(dVal)
        match(k) = true;
    end
end

filesD  = files_e(match);
aValsD  = aVals(match);
bValsD  = bVals(match);
[aValsD, idx] = sort(aValsD);
filesD  = filesD(idx);
bValsD  = bValsD(idx);

keep   = arrayfun(@(a) any(abs(a - [2,3,4]) < 1e-9), aValsD);
filesD = filesD(keep); aValsD = aValsD(keep); bValsD = bValsD(keep);

if isempty(filesD)
    warning('Panel (b): no files found for a=2,3,4 at d=%.6f.', targetD_ellipse);
else
    sepHandles   = gobjects(numel(filesD),1);
    legendEntries = strings(numel(filesD),1);
    xStagger = targetD_ellipse + 3.0*max(aValsD);

    for k = 1:numel(filesD)
        filePath = fullfile(filesD(k).folder, filesD(k).name);
        A = readmatrix(filePath);
        if isempty(A) || size(A,2) < 8; continue; end

        t  = A(:,1);
        x1 = A(:,2); y1 = A(:,3); th1 = A(:,4);
        x2 = A(:,5); y2 = A(:,6); th2 = A(:,7);
        bbDist = A(:,8);

        a = aValsD(k); b = bValsD(k);
        cLine = cmap(round(a),:);  % a=2,3,4 → row 2,3,4
        sepHandles(k) = plot(axSep, t, bbDist, 'Color', cLine);

        legendEntries(k) = strcat('{\it a/b} = ', num2str(round(a/b)));

        nSnap = min(10, numel(t));
        tSnap = linspace(t(1), t(end), nSnap);
        x1s = interp1(t, x1, tSnap); y1s = interp1(t, y1, tSnap);
        x2s = interp1(t, x2, tSnap); y2s = interp1(t, y2, tSnap);
        th1s = interp1(t, unwrap(th1), tSnap);
        th2s = interp1(t, unwrap(th2), tSnap);
        xShift = (k-1)*xStagger;
        phi_e = linspace(0, 2*pi, 80);
        ex0 = a*cos(phi_e); ey0 = b*sin(phi_e);

        for j = 1:nSnap
            alpha = 0.15 + 0.85*(j-1)/max(1, nSnap-1);
            cFill = blend_with_white(blend_with_white(cLine, alpha), 0.9);
            ct1 = cos(th1s(j)); st1 = sin(th1s(j));
            patch(axInsetSep, x1s(j)+ct1*ex0-st1*ey0+xShift, y1s(j)+st1*ex0+ct1*ey0, ...
                cFill, 'FaceColor', cFill, 'EdgeColor', 'k', 'LineWidth', 0.7);
            ct2 = cos(th2s(j)); st2 = sin(th2s(j));
            patch(axInsetSep, x2s(j)+ct2*ex0-st2*ey0+xShift, y2s(j)+st2*ex0+ct2*ey0, ...
                cFill, 'FaceColor', cFill, 'EdgeColor', 'k', 'LineWidth', 0.7);
        end
    end

    xlabel(axSep, 'time, {\it t}');
    ylabel(axSep, 'dist({\it \Gamma}_1, {\it \Gamma}_2)');
    valid = isgraphics(sepHandles) & strlength(legendEntries) > 0;
    legend(axSep, sepHandles(valid), cellstr(legendEntries(valid)), 'Location', 'northeast', 'Color', 'none');
    axis(axInsetSep, 'equal');
    axInsetSep.XTick = []; axInsetSep.YTick = [];
end

% ==========================================================================
% Layout
% ==========================================================================
axCirc.ActivePositionProperty = 'position';
axSep.ActivePositionProperty  = 'position';
axCirc.Units = 'inches'; axSep.Units = 'inches';
axInsetCirc.Units = 'inches'; axInsetSep.Units = 'inches';

tiC = axCirc.TightInset; tiS = axSep.TightInset;

leftMargin   = max([0.8, tiC(1), tiS(1)]);
rightMargin  = max([0.22, tiC(3), tiS(3)]) + 0.06;
bottomMargin = max([0.42, tiC(2), tiS(2)]) + 0.08;
topMargin    = max([0.24, tiC(4), tiS(4)]) + 0.06;
midGap       = 1.2;

panelW = (fig.PaperSize(1) - leftMargin - rightMargin - midGap) / 2;
panelH = fig.PaperSize(2) - bottomMargin - topMargin;

axCirc.Position = [leftMargin,                   bottomMargin, panelW, panelH];
axSep.Position  = [leftMargin + panelW + midGap, bottomMargin, panelW, panelH];

% Insets: fill most of each panel, anchored to panel lower-left.
axInsetCirc.ActivePositionProperty = 'position';
axInsetCirc.Position = [axCirc.Position(1) + 0.45, axCirc.Position(2) + 0.8, 0.4*panelW, 0.4*panelH];
axInsetSep.Position  = [axSep.Position(1)  - 0.7,         axSep.Position(2)  - 0.7,          0.8*panelW, 0.8*panelH];

% Panel labels.
a_ann = panel_label('\it (a)', axCirc);
b_ann = panel_label('\it (b)', axSep);

% ---- helpers ---------------------------------------------------------------
function val = parse_param(fname, key)
    token = regexp(fname, ['_' key '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'], 'tokens', 'once');
    if isempty(token); val = NaN; else; val = str2double(token{1}); end
end

function c = blend_with_white(baseColor, alphaLike)
    alphaLike = max(0, min(1, alphaLike));
    c = alphaLike * baseColor + (1 - alphaLike) * [1, 1, 1];
end

function a = panel_label(str, ax)
    a = annotation('textbox', 'String', str, 'FontName', 'Times', ...
        'EdgeColor', 'none', 'FontSize', 18, 'Units', 'inches');
    a.Position(1) = ax.Position(1) - 0.8;
    a.Position(2) = sum(ax.Position([2 4])) - a.Position(4) + 0.2;
end