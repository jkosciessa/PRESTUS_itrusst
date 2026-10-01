function moved = itrusst_backup_outputs(rootpath, tag)
% Move all simulation outputs and figures aside before a rerun.
%
% Keeps the I/O structure: after the move, data/tussim/, data/calibration/
% and figures/ are empty and the scripts write to the same paths as before.
%
%   data/tussim       -> data/tussim_<tag>
%   data/calibration  -> data/calibration_<tag>          (if present)
%   figures/*         -> figures/backup/figures_<tag>/*  (except figures/backup)
%
% Moves are renames on the same volume (no copy). Refuses to overwrite an
% existing backup with the same tag.
%
% Inputs:
%   rootpath - project root
%   tag      - backup suffix, e.g. '261001' (yymmdd, as tussim_260520)
% Output:
%   moved    - n x 2 cell array {source, target} of the moves made

pn.data    = fullfile(rootpath, 'data');
pn.figures = fullfile(rootpath, 'figures');

pairs = {fullfile(pn.data, 'tussim'),      fullfile(pn.data, ['tussim_' tag]); ...
         fullfile(pn.data, 'calibration'), fullfile(pn.data, ['calibration_' tag])};

fig_backup = fullfile(pn.figures, 'backup', ['figures_' tag]);
fig_items = dir(pn.figures);
fig_items = fig_items(~ismember({fig_items.name}, {'.', '..', 'backup', '.DS_Store'}));
for i = 1:numel(fig_items)
    pairs(end+1, :) = {fullfile(pn.figures, fig_items(i).name), fullfile(fig_backup, fig_items(i).name)}; %#ok<AGROW>
end

% check all targets first, so that a clash does not leave a half-moved state
pairs = pairs(cellfun(@(p) exist(p, 'file') > 0, pairs(:, 1)), :);
clash = cellfun(@(p) exist(p, 'file') > 0, pairs(:, 2));
if any(clash)
    error('Backup target already exists: %s', strjoin(pairs(clash, 2), ', '));
end

if ~isempty(fig_items) && ~isfolder(fig_backup); mkdir(fig_backup); end
for i = 1:size(pairs, 1)
    [ok, msg] = movefile(pairs{i, 1}, pairs{i, 2});
    if ~ok
        error('Could not move %s -> %s: %s', pairs{i, 1}, pairs{i, 2}, msg);
    end
    fprintf('moved %s -> %s\n', pairs{i, 1}, pairs{i, 2});
end
mkdir(fullfile(pn.data, 'tussim'));

moved = pairs;
