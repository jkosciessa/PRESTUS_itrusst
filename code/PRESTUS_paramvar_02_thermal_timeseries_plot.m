% PRESTUS_PARAMVAR_02_THERMAL_TIMESERIES_PLOT
%
% Plots the full thermal time-course (not just the peak/end-of-sonication
% scalars used elsewhere) at two locations, overlaid in one plot:
%
%   Focus      - temperature at the empirical acoustic focus voxel
%                (Ix_brain_vox, Iz_brain_vox; per condition, see
%                load_focus_heating_metrics.m)
%   Peak skull - at each time point, the hottest skull voxel (skull mask
%                from final_tissues.nii.gz, label==4); re-evaluated at
%                every time step, not fixed to one location
%
% Both locations are overlaid in a single plot: for each, the min-max
% range across ALL paramvar sweep conditions (28 tissue x property x
% {bottom,top} combinations, plus the mid-level baseline) is drawn as a
% shaded band, with the across-condition mean overlaid - the
% time-resolved analogue of the scalar range plots in
% plot_metric_with_patches.m / PRESTUS_paramvar_02_benchmark2_plot.m.
%
% Time-resolved data come from the thermal cache .mat files
% (sub-XXX_phantom_heating_res<affix>.mat -> results_heating.focal_planeT,
% time_status_seq), written when parameters.io.save_thermal_matrices is
% enabled (the default). Reading + reducing all 29 of these per run is
% slow, so the reduced [n_cond x n_time] matrices are cached to a single
% aggregate .mat (agg_file below) after the first run; set force_reload
% to true to rebuild it from the raw cache files.
%
% Output: figures/paramvar/sub-XXX_thermal_timeseries.png/.eps

%% Paths

currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code'));

pn.figures = fullfile(rootpath, 'figures', 'paramvar');
if ~isfolder(pn.figures); mkdir(pn.figures); end
pn.simnibs = fullfile(rootpath, 'data', 'simnibs');

%% Settings

subject_id   = 002;
variant_name = 'itrusst_protocol1a';
sub_str      = sprintf('sub-%03d', subject_id);
pn.data      = fullfile(rootpath, 'data', 'tussim', variant_name, sub_str);

tissue_labels      = {'brain', 'skull'};
tissue_suffixes    = {'_b', '_s'};
property_labels    = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                      'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes  = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};
magnitude_suffixes = {'1', '2'}; % bottom, top

force_reload = false; % set true to ignore the cached aggregate and reload from the raw per-condition cache files
agg_file = fullfile(pn.data, 'cache', sprintf('%s_phantom_thermal_timeseries_agg.mat', sub_str));

% Build the list of all condition affixes: mid-level baseline + the 28
% bottom/top sweep conditions.
affixes = {'_baseline'};
labels  = {'baseline'};
for t = 1:numel(tissue_labels)
    for p = 1:numel(property_labels)
        for v = 1:numel(magnitude_suffixes)
            affixes{end+1} = [tissue_suffixes{t} property_suffixes{p} magnitude_suffixes{v}]; %#ok<SAGROW>
            labels{end+1}  = [tissue_labels{t}(1) '_' property_labels{p} '_' magnitude_suffixes{v}]; %#ok<SAGROW>
        end
    end
end

%% Load per-condition time series (or reload the cached aggregate)

if isfile(agg_file) && ~force_reload
    fprintf('Loading cached aggregate: %s\n', agg_file);
    load(agg_file, 't_ref', 'focus_ts', 'skull_ts', 'loaded_labels');
else
    % Skull mask, shared across all conditions (fixed phantom geometry;
    % only tissue properties vary across the sweep, not tissue placement)
    seg_path = fullfile(pn.simnibs, sprintf('m2m_%s', sub_str), 'final_tissues.nii.gz');
    seg      = round(double(niftiread(seg_path)));   % [Nx x Nz], labels: 0 water, 1 brain, 4 skull
    skull_mask = (seg == 4);
    if ~any(skull_mask(:))
        warning('No skull tissue found in %s - "peak skull" trace will be empty.', seg_path);
    end

    t_ref        = [];   % common reference time axis, taken from the first loaded condition
    focus_ts     = [];   % [n_cond x n_time], deg C rise from 37
    skull_ts     = [];
    loaded_labels = {};  % condition label for each row of focus_ts/skull_ts

    for c = 1:numel(affixes)
        affix    = affixes{c};
        csv_file = fullfile(pn.data, [sub_str '_phantom' affix '.csv']);
        mat_file = fullfile(pn.data, 'cache', [sub_str '_phantom_heating_res' affix '.mat']);

        if ~isfile(csv_file) || ~isfile(mat_file)
            warning('Skipping %s - missing csv or cache file.', labels{c});
            continue;
        end

        try
            Tcsv     = readtable(csv_file);
            ix_focus = Tcsv.Ix_brain_vox(1);
            iz_focus = Tcsv.Iz_brain_vox(1);

            S        = load(mat_file, 'results_heating', 'time_status_seq');
            recorded = logical([S.time_status_seq.recorded]);
            t_sec    = [S.time_status_seq.time];
            t_sec    = t_sec(recorded);

            planeT = S.results_heating.focal_planeT(:, :, recorded);   % [Nx x Nz x n_time]
            [Nx, Nz, nt] = size(planeT);

            this_focus = squeeze(planeT(ix_focus, iz_focus, :))' - 37;

            if any(skull_mask(:))
                flatT      = reshape(planeT, Nx * Nz, nt);
                this_skull = max(flatT(skull_mask(:), :), [], 1) - 37;
            else
                this_skull = nan(1, nt);
            end
        catch err
            warning('Skipping %s - %s', labels{c}, err.message);
            continue;
        end

        if isempty(t_ref)
            t_ref = t_sec;
        elseif numel(t_sec) ~= numel(t_ref) || any(abs(t_sec - t_ref) > 1e-6)
            % Interpolate onto the reference time axis if a condition's time
            % vector differs (should not normally happen - same protocol).
            this_focus = interp1(t_sec, this_focus, t_ref, 'linear', 'extrap');
            this_skull = interp1(t_sec, this_skull, t_ref, 'linear', 'extrap');
        end

        focus_ts(end+1, :)      = this_focus; %#ok<AGROW>
        skull_ts(end+1, :)      = this_skull; %#ok<AGROW>
        loaded_labels{end+1}    = labels{c};   %#ok<AGROW>
    end

    fprintf('Loaded %d/%d conditions.\n', size(focus_ts, 1), numel(affixes));

    if ~isfolder(fileparts(agg_file)); mkdir(fileparts(agg_file)); end
    save(agg_file, 't_ref', 'focus_ts', 'skull_ts', 'loaded_labels', '-v7.3');
    fprintf('Saved aggregate to %s\n', agg_file);
end

%% Plot: both locations overlaid, mean +- range across all paramvar settings

TRACES = struct( ...
    'data',  {focus_ts,          skull_ts}, ...
    'name',  {'Focus (brain)',   'Peak skull'}, ...
    'color', {[0.35 0.35 0.35],  [0.7 0 0]});

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 20 9]);
set(fig, 'Renderer', 'painters');
ax = axes('Parent', fig); hold(ax, 'on');

% Pass 1: shaded min-max bands (drawn first, so mean lines stay on top)
for i_trace = 1:numel(TRACES)
    data = TRACES(i_trace).data;
    if isempty(data) || all(isnan(data(:))); continue; end
    clr      = TRACES(i_trace).color;
    band_clr = 1 - (1 - clr) * 0.20;
    lo = min(data, [], 1);
    hi = max(data, [], 1);
    fill(ax, [t_ref, fliplr(t_ref)], [lo, fliplr(hi)], band_clr, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Pass 2: across-condition mean lines
for i_trace = 1:numel(TRACES)
    data = TRACES(i_trace).data;
    if isempty(data) || all(isnan(data(:))); continue; end
    clr = TRACES(i_trace).color;
    mn  = mean(data, 1, 'omitnan');
    plot(ax, t_ref, mn, '-', 'Color', clr, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('%s (n=%d)', TRACES(i_trace).name, size(data, 1)));
end

xlabel(ax, 'Time (s)');
ylabel(ax, '\DeltaT (\circC)', 'Interpreter', 'tex');
title(ax, sprintf('%s - thermal time-course: mean \\pm range across paramvar settings', sub_str), ...
    'Interpreter', 'tex');
legend(ax, 'show', 'Location', 'northwest', 'Box', 'off');
box(ax, 'on');

base = fullfile(pn.figures, sprintf('%s_thermal_timeseries', sub_str));
print(fig, [base '.png'], '-dpng', '-r300');
saveas(fig, [base '.eps'], 'epsc');
close(fig);
fprintf('Saved %s\n', base);
