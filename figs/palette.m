function P = palette()
% PALETTE  Return a struct of colors and a colormap sampled from a cmocean map.
%
% In dark mode (makefigs dark), foreground and background colors are swapped and curves are
% colored from the lighter part of the map, so they remain visible on a dark page.
name   = 'dense';
raw    = cmocean(name, 7);
P.cmap = @(varargin) cmocean(name, varargin{:});   % P.cmap() or P.cmap(n)
P.c1   = raw(7,:);
P.c2   = raw(3,:);
P.col  = raw(5,:);
P.lines = raw(2:5,:);
P.snap  = raw(3:5,:);

P.dark    = isappdata(groot, 'figdark') && getappdata(groot, 'figdark');
P.linemap = P.cmap;                % colors for curves, P.linemap(n)
P.fg      = [0 0 0];               % text, axes and lines on the page
P.bg      = [1 1 1];               % page, e.g. marker faces
P.ink     = [0.2 0.2 0.2];         % secondary lines
P.grey    = [0.4 0.4 0.4];         % trajectories
P.shade   = [0.9 0.9 0.9];         % shaded regions

if P.dark
    full      = cmocean(name, 256);
    P.linemap = @(n) interp1(linspace(0, 1, 256), full, linspace(0, 0.62, n));
    P.c1      = P.linemap(1);
    P.fg      = [0.92 0.90 0.94];
    P.bg      = [23 18 27] / 255;  % the presentation's dark page colour
    P.ink     = [0.75 0.73 0.78];
    P.grey    = [0.6 0.58 0.63];
    P.shade   = [0.25 0.22 0.28];
end
end
