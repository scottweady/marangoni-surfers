
% makefigs.m  —  Run figure scripts and save PDFs to figs/output/
%
%   makefigs                   all figures
%   makefigs single pair       only the named groups: single, pair, suspension
%   makefigs fig3 fig8         only the named figures, numbered as in the manuscript
function makefigs(varargin)

figsDir = fileparts(mfilename('fullpath'));
addpath(figsDir);
dataRoot = fullfile(figsDir, '..', 'data');
outDir  = fullfile(figsDir, 'output');
if ~exist(outDir, 'dir'), mkdir(outDir); end

groups = struct( ...
    'single',     {{'fig_schematic', 'fig_multipole', 'fig_convergence'}}, ...
    'pair',       {{'fig_pair_parallel', 'fig_pair_antiparallel', 'fig_equilibrium_separation'}}, ...
    'suspension', {{'fig_expansion', 'fig_patches_farfield', 'fig_tandem'}});

% figures in the order they appear in the manuscript, selected as fig1, fig2, ...
numbered = {'fig_schematic', 'fig_multipole', 'fig_pair_parallel', 'fig_pair_antiparallel', ...
            'fig_equilibrium_separation', 'fig_expansion', 'fig_patches_farfield', ...
            'fig_tandem', 'fig_convergence'};

selected = varargin;
if isempty(selected), selected = fieldnames(groups)'; end
scripts = cell(1, numel(selected));
for k = 1:numel(selected)
    n = sscanf(selected{k}, 'fig%d');
    if isfield(groups, selected{k})
        scripts{k} = groups.(selected{k});
    elseif ~isempty(n) && strcmp(selected{k}, sprintf('fig%d', n)) && n >= 1 && n <= numel(numbered)
        scripts{k} = numbered(n);
    else
        error('Unknown figure: %s (expected single, pair, suspension, or fig1 to fig%d)', ...
              selected{k}, numel(numbered));
    end
end
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
