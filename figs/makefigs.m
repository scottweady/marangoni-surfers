
% makefigs.m  —  Run figure scripts and save PDFs to figs/output/
%
%   makefigs                   all figures
%   makefigs single pair       only the named groups: single, pair, suspension
function makefigs(varargin)

figsDir = fileparts(mfilename('fullpath'));
addpath(figsDir);
dataRoot = fullfile(figsDir, '..', 'data');
outDir  = fullfile(figsDir, 'output');
if ~exist(outDir, 'dir'), mkdir(outDir); end

groups = struct( ...
    'single',     {{'fig_schematic', 'fig_multipole', 'fig_stability', 'fig_convergence'}}, ...
    'pair',       {{'fig_pair_parallel', 'fig_pair_antiparallel', 'fig_equilibrium_separation'}}, ...
    'suspension', {{'fig_expansion', 'fig_patches'}});

selected = varargin;
if isempty(selected), selected = fieldnames(groups)'; end
unknown = setdiff(selected, fieldnames(groups));
if ~isempty(unknown)
    error('Unknown figure group: %s (expected single, pair, or suspension)', strjoin(unknown, ', '));
end
scripts = cellfun(@(g) groups.(g), selected, 'UniformOutput', false);
scripts = [scripts{:}];

% A figure whose script fails (e.g. its data has not been generated) is skipped
ok = cellfun(@(name) make_figure(name, figsDir, dataRoot, outDir), scripts);

fprintf('Done. PDFs written to %s\n', outDir);
if any(~ok)
    fprintf('Skipped: %s\n', strjoin(scripts(~ok), ', '));
end

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
