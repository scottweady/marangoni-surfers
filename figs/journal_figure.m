
function fig = journal_figure(sz, sc, tex)
% JOURNAL_FIGURE  Create a figure formatted for journal publication.

% Set scale
if nargin < 2
  sc = 1;
end

% Set interpreter
if nargin < 3
  tex = 'latex';
end

% Initialize
fig = figure;
fig.Units = 'inches';
fig.PaperSize = sz * sc;
fig.PaperPosition = [0 0 sz] * sc;
fig.Position = [2 2 sz * sc];
fig.PaperPositionMode = 'auto';

% Set defaults
set(fig,'defaultlinelinewidth', 1.25 * sc)
set(fig,'defaultAxesLineWidth', 1 * sc)
set(fig,'defaultAxesFontSize', 10 * sc)
set(fig,'defaulttextinterpreter',tex); 
set(fig,'defaultAxesTickLabelInterpreter',tex); 
set(fig,'defaultLegendInterpreter',tex);
set(fig,'defaultLegendFontSize', 10 * sc);
set(fig, 'defaultLegendEdgeColor', 'none');
set(fig,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'});
set(fig,'DefaultTextColor','k');
set(fig,'defaultAxesFontName','Times');
set(fig,'defaultTextFontName','Times');

% Dark mode (makefigs dark): light text and axes on the dark page colour. The background is
% opaque, since MATLAB paints a transparent one white wherever a figure uses transparency.
P = palette();
if P.dark
  fig.Color = P.bg;
  fig.InvertHardcopy = 'off';
  set(fig,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{P.fg,P.fg,P.fg});
  set(fig,'DefaultAxesColor',P.bg);
  set(fig,'DefaultTextColor',P.fg);
  set(fig,'DefaultLegendTextColor',P.fg);
  set(fig,'DefaultColorbarColor',P.fg);
  set(fig,'DefaultTextboxshapeColor',P.fg);
end

