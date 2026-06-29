
% make_figs.m  —  Run all figure scripts and save PDFs to figs/output/
clear

figsDir = fileparts(mfilename('fullpath'));
addpath(figsDir);
P = palette();   % shared color palette (cmocean 'tempo')
outDir  = fullfile(figsDir, 'output');
if ~exist(outDir, 'dir'), mkdir(outDir); end

scripts = { 'fig_pair_prolate', 'fig_pair_disk', 'fig_multipole', 'fig_center_speed', 'fig_expansion', 'fig_pair_oblate' };

for k = 1:numel(scripts)
    name = scripts{k};
    fprintf('Running %s ...\n', name);
    run(fullfile(figsDir, name));
    print(gcf, fullfile(outDir, name), '-dpdf', '-vector');
    fprintf('  -> %s.pdf\n', name);
    close all;
end

fprintf('Done. PDFs written to %s\n', outDir);
