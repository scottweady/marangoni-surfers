% thin_lines.m — keep every n-th vertex of long lines (e.g. streamslice output, which is
% sampled far finer than print resolution), preserving endpoints and NaN breaks.
function thin_lines(h, n)
for hk = reshape(h, 1, [])
    x = hk.XData;  y = hk.YData;
    m = numel(x);
    if m <= 10*n, continue; end
    keep = false(1, m);
    keep([1:n:m, m]) = true;
    gap = find(isnan(x));
    keep([gap, max(gap - 1, 1), min(gap + 1, m)]) = true;
    set(hk, 'XData', x(keep), 'YData', y(keep));
end
end
