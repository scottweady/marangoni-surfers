function P = palette()
% TEMPO_PALETTE  Return a struct of colors sampled from cmocean('tempo').
raw    = cmocean('tempo', 7);
P.c1   = raw(7,:);
P.c2   = raw(3,:);
P.col  = raw(5,:);
P.lines = raw(2:5,:);
P.snap  = raw(3:5,:);
end
