
% make_figs.m  —  Run all figure scripts and save PDFs to figs/output/
clear

figsDir = fileparts(mfilename('fullpath'));
addpath(figsDir);
dataRoot = fullfile(figsDir, '..', 'data');
P = palette();   % shared color palette (cmocean 'tempo')
outDir  = fullfile(figsDir, 'output');
if ~exist(outDir, 'dir'), mkdir(outDir); end

scripts = {'fig_schematic', ...
'fig_multipole', ...
'fig_pair_parallel', ...
'fig_pair_antiparallel', ...
'fig_equilibrium_separation', ...
'fig_expansion', ...
'fig_patches', ...
'fig_stability', ...
'fig_convergence'};

for k = 1:numel(scripts)
    name = scripts{k};
    fprintf('Running %s ...\n', name);
    run(fullfile(figsDir, name));
    print(gcf, fullfile(outDir, name), '-dpdf', '-vector');
    fprintf('  -> %s.pdf\n', name);
    close all;
end

fprintf('Done. PDFs written to %s\n', outDir);
