
% make_figs.m  —  Run all figure scripts and save PDFs to figs/output/
clear

figsDir = fileparts(mfilename('fullpath'));
addpath(figsDir);
dataRoot = fullfile(figsDir, '..', 'data-07-10-2026');
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

% A figure whose script fails (e.g. its data has not been generated) is skipped
ok = cellfun(@(name) make_figure(name, figsDir, dataRoot, outDir), scripts);

fprintf('Done. PDFs written to %s\n', outDir);
if any(~ok)
    fprintf('Skipped: %s\n', strjoin(scripts(~ok), ', '));
end

% Run one figure script in its own workspace, so its variables cannot overwrite the
% driver's, and print it to PDF. Returns false if the script fails.
function ok = make_figure(name, figsDir, dataRoot, outDir)   %#ok<INUSL> dataRoot is read by the scripts
    P = palette();   %#ok<NASGU> shared color palette (cmocean 'tempo'), used by the scripts
    pdfFile = fullfile(outDir, name);
    fprintf('Running %s ...\n', name);
    try
        run(fullfile(figsDir, name));
        print(gcf, pdfFile, '-dpdf', '-vector');
        fprintf('  -> %s.pdf\n', pdfFile);
        ok = true;
    catch err
        fprintf('  skipped: %s\n', err.message);
        ok = false;
    end
    close all;
end
