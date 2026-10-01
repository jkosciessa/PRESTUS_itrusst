% CALIBRATION_BENCHMARK_02_PLOT
% Compare calibration outputs between the 3D and axisymmetric (2D) variants
% produced by explore_calibration_benchmark_01_sim.m.
%
% Reads the YAML files written by save_optimized_values and the benchmark_results.mat
% timing summary, then produces:
%   (1) Bar chart of wall time per variant
%   (2) Comparison of optimised source amplitude per intensity
%   (3) Comparison of optimised element phases
%   (4) Difference in amplitude and phase between variants

clear all; close all; clc;

% get root path
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, 'julkos');

%% Configuration — must match explore_calibration_benchmark_01_sim.m

variant_names   = {'3D', 'axisymmetric2D'};
variant_labels  = {'3D', '2D axisymmetric'};
bench_dir       = fullfile(rootpath, 'data', 'calibration_benchmark');
fig_dir         = fullfile(rootpath, 'figures', 'debug', 'calibration_benchmark');
if ~isfolder(fig_dir); mkdir(fig_dir); end

%% Load timing summary

timing_file = fullfile(bench_dir, 'benchmark_results.mat');
if ~isfile(timing_file)
    error('Run explore_calibration_benchmark_01_sim.m first to generate %s', timing_file);
end
load(timing_file, 'results');

%% Load calibrated YAML files for each variant

% Discover YAML files from the first variant (both should share the same set)
profile_dir_ref = fullfile(bench_dir, variant_names{1}, 'PRESTUS_virtual_parameters');
yaml_files = dir(fullfile(profile_dir_ref, '*.yaml'));

if isempty(yaml_files)
    error('No calibrated YAML files found in %s', profile_dir_ref);
end

n_files    = numel(yaml_files);
n_variants = numel(variant_names);

% Pre-allocate storage indexed by file
cal = struct();
cal.yaml_names    = {yaml_files.name};
cal.focal_dist_mm = nan(1, n_files);
cal.intensity     = nan(1, n_files);
cal.elem_amp      = nan(n_variants, n_files);
cal.elem_phase_deg = cell(n_variants, n_files);  % cell because n_elements may vary

for f = 1:n_files
    fname = yaml_files(f).name;

    % Parse focal distance and intensity from filename:
    %   {equip}-F{focal}mm-I{intens}wpercm2.yaml
    tok_F = regexp(fname, '-F([\d.]+)mm', 'tokens', 'once');
    tok_I = regexp(fname, '-I([\d.]+)wpercm2', 'tokens', 'once');
    if isempty(tok_F) || isempty(tok_I)
        warning('Cannot parse focal/intensity from filename: %s — skipping.', fname);
        continue
    end
    cal.focal_dist_mm(f) = str2double(tok_F{1});
    cal.intensity(f)     = str2double(tok_I{1});

    for v = 1:n_variants
        yaml_path = fullfile(bench_dir, variant_names{v}, 'PRESTUS_virtual_parameters', fname);
        if ~isfile(yaml_path)
            warning('Missing YAML for variant %s: %s', variant_names{v}, fname);
            continue
        end
        data = yaml.loadFile(yaml_path);
        cal.elem_amp(v, f)       = data.transducer.annular.elem_amp(1);
        cal.elem_phase_deg{v, f} = data.transducer.annular.elem_phase_deg;
    end
end

% Sort files by intensity then focal distance for consistent plotting
[~, sort_idx] = sortrows([cal.intensity(:), cal.focal_dist_mm(:)]);
cal.focal_dist_mm  = cal.focal_dist_mm(sort_idx);
cal.intensity      = cal.intensity(sort_idx);
cal.elem_amp       = cal.elem_amp(:, sort_idx);
cal.elem_phase_deg = cal.elem_phase_deg(:, sort_idx);
cal.yaml_names     = cal.yaml_names(sort_idx);
n_files = numel(sort_idx);

%% x-axis labels

xlabels = arrayfun(@(f) sprintf('F%.0f mm\nI=%.0f W/cm²', ...
    cal.focal_dist_mm(f), cal.intensity(f)), 1:n_files, 'UniformOutput', false);

colors = lines(n_variants);

%% Figure 1: Wall time

h1 = figure('Name', 'Calibration wall time', 'Position', [50 50 500 350]);
wall_times = [results.wall_time_s];
bar(wall_times / 60, 0.5, 'FaceColor', [0.4 0.6 0.9]);
xticks(1:n_variants);
xticklabels(variant_labels);
ylabel('Wall time (min)');
title('Calibration wall time');
for v = 1:n_variants
    text(v, wall_times(v)/60 + 0.02*max(wall_times/60), ...
        sprintf('%.1f min', wall_times(v)/60), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end
saveas(h1, fullfile(fig_dir, 'walltime_comparison.png'), 'png');
saveas(h1, fullfile(fig_dir, 'walltime_comparison.eps'), 'epsc');

%% Figure 2: Source amplitude per condition

h2 = figure('Name', 'Source amplitude', 'Position', [50 420 700 350]);
hold on
for v = 1:n_variants
    plot(1:n_files, cal.elem_amp(v,:), 'o-', ...
        'Color', colors(v,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', variant_labels{v});
end
xticks(1:n_files);
xticklabels(xlabels);
ylabel('Source amplitude (Pa)');
title('Optimised source amplitude');
legend('Location', 'best');
grid on;
saveas(h2, fullfile(fig_dir, 'amplitude_comparison.png'), 'png');
saveas(h2, fullfile(fig_dir, 'amplitude_comparison.eps'), 'epsc');

%% Figure 3: Amplitude difference (3D as reference)

h3 = figure('Name', 'Amplitude difference', 'Position', [760 420 700 350]);
amp_diff = cal.elem_amp(2,:) - cal.elem_amp(1,:);
amp_diff_pct = 100 * amp_diff ./ cal.elem_amp(1,:);
bar(amp_diff_pct, 0.5, 'FaceColor', [0.8 0.4 0.3]);
xticks(1:n_files);
xticklabels(xlabels);
yline(0, 'k--');
ylabel('Amplitude difference (%)');
title('2D axisymmetric minus 3D (% of 3D)');
grid on;
saveas(h3, fullfile(fig_dir, 'amplitude_difference.png'), 'png');
saveas(h3, fullfile(fig_dir, 'amplitude_difference.eps'), 'epsc');

%% Figure 4: Element phases (one sub-panel per condition)

n_elems = max(cellfun(@numel, cal.elem_phase_deg(1,:)));
if n_elems > 1
    h4 = figure('Name', 'Element phases', 'Position', [50 800 900 300*ceil(n_files/3)]);
    for f = 1:n_files
        subplot(ceil(n_files/3), 3, f);
        hold on;
        for v = 1:n_variants
            ph = cal.elem_phase_deg{v,f};
            if ~isempty(ph)
                plot(1:numel(ph), ph, 'o-', ...
                    'Color', colors(v,:), 'DisplayName', variant_labels{v});
            end
        end
        xticks(1:n_elems);
        xlabel('Element');
        ylabel('Phase (°)');
        title(xlabels{f}, 'Interpreter', 'none');
        if f == 1; legend('Location', 'best'); end
        grid on;
    end
    saveas(h4, fullfile(fig_dir, 'phase_comparison.png'), 'png');
    saveas(h4, fullfile(fig_dir, 'phase_comparison.eps'), 'epsc');
end

%% Console summary

fprintf('\n%-20s  %10s  %10s  %10s\n', 'Condition', '3D amp', '2D amp', 'Diff (%)');
fprintf('%s\n', repmat('-', 1, 55));
for f = 1:n_files
    fprintf('%-20s  %10.0f  %10.0f  %+9.1f%%\n', ...
        xlabels{f}, cal.elem_amp(1,f), cal.elem_amp(2,f), amp_diff_pct(f));
end

fprintf('\n%-20s  %8s\n', 'Variant', 'Time (min)');
fprintf('%s\n', repmat('-', 1, 30));
for v = 1:n_variants
    fprintf('%-20s  %8.1f\n', variant_labels{v}, results(v).wall_time_s/60);
end

fprintf('\nFigures saved to %s\n', fig_dir);
