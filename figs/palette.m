function P = palette()
% PALETTE  Return a struct of colors and a colormap sampled from a cmocean map.
name   = 'dense';
raw    = cmocean(name, 7);
P.cmap = @(varargin) cmocean(name, varargin{:});   % P.cmap() or P.cmap(n)
P.c1   = raw(7,:);
P.c2   = raw(3,:);
P.col  = raw(5,:);
P.lines = raw(2:5,:);
P.snap  = raw(3:5,:);
end
