% Extract post-skull voxel metrics for sub-002.
% Reads Isppa at a fixed post-skull voxel (35,35,58) from the intensity NIfTI
% and riseT_brain from the output-table CSV (no heatrise NIfTI is produced).

% get root path (script must be run from file location)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.figures  = fullfile(rootpath, 'figures');
pn.data     = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a');
pn.code     = fullfile(rootpath, 'code'); addpath(pn.code);

subdir_path  = fullfile(pn.data, 'sub-002');
subject_name = 'sub-002';
postskull_vox = [35, 35, 58]; % fixed post-skull voxel (x,y,z); 2D NIfTI indexed as V(x,z)

% Define tissues and suffixes
tissue_labels   = {'brain', 'skull'};
tissue_suffixes = {'_b', '_s'};

% Property definitions (must match simulation script)
property_labels   = {'sound_speed', 'density', 'alpha_coeff', ...
                     'thermal_conductivity', 'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};
property_symbols  = {'c', '\rho', '\alpha', '\kappa', 'c_p', '\omega', '\beta'};

sweep_values   = 1:2;
scatter_levels = [1, 2];
scatter_colors = {[1 1 1], [0 0 0]};

n_conds = length(tissue_labels) * length(property_labels);
isppa_postskull    = nan(length(sweep_values), n_conds);
heatrise_postskull = nan(length(sweep_values), n_conds);
xlabels = {};

i_cond = 0;
for t = 1:length(tissue_labels)
    tissue_suf = tissue_suffixes{t};

    for p = 1:length(property_labels)
        i_cond = i_cond + 1;
        prop_label = property_labels{p};
        prop_suf   = property_suffixes{p};
        xlabels{i_cond} = [property_symbols{p} '_{' tissue_suf(2:end) '}'];

        for i_sw = 1:length(sweep_values)
            sw       = sweep_values(i_sw);
            affix    = [tissue_suf prop_suf num2str(sw)];
            csv_file = fullfile(subdir_path, ['sub-002_phantom' affix '.csv']);

            if ~isfile(csv_file)
                warning('CSV not found: %s', csv_file);
                continue
            end

            T_csv = readtable(csv_file);

            % Isppa at the fixed post-skull voxel from the intensity NIfTI
            intensity_pat   = fullfile(subdir_path, 'nii', sprintf('*phantom_T1w%s%s%d_intensity.nii.gz', tissue_suf, prop_suf, sw));
            intensity_files = dir(intensity_pat);
            if ~isempty(intensity_files)
                try
                    V = niftiread(fullfile(intensity_files(1).folder, intensity_files(1).name));
                    isppa_postskull(i_sw, i_cond) = V(postskull_vox(1), postskull_vox(3));
                    fprintf('%s sw%d (%s %s): Isppa voxel(%d,%d)=%.4f\n', ...
                        subject_name, sw, tissue_labels{t}, prop_label, ...
                        postskull_vox(1), postskull_vox(3), V(postskull_vox(1), postskull_vox(3)));
                catch ME
                    warning('Intensity NIfTI load failed %s: %s', intensity_files(1).name, ME.message);
                end
            else
                warning('Intensity NIfTI not found: %s', intensity_pat);
            end

            % Heatrise at the fixed post-skull voxel from the heatrise NIfTI
            heatrise_pat   = fullfile(subdir_path, 'nii', sprintf('*phantom_T1w%s%s%d_heatrise.nii.gz', tissue_suf, prop_suf, sw));
            heatrise_files = dir(heatrise_pat);
            if ~isempty(heatrise_files)
                try
                    Vh = niftiread(fullfile(heatrise_files(1).folder, heatrise_files(1).name));
                    heatrise_postskull(i_sw, i_cond) = Vh(postskull_vox(1), postskull_vox(3));
                    fprintf('%s sw%d (%s %s): heatrise voxel(%d,%d)=%.4f\n', ...
                        subject_name, sw, tissue_labels{t}, prop_label, ...
                        postskull_vox(1), postskull_vox(3), Vh(postskull_vox(1), postskull_vox(3)));
                catch ME
                    warning('Heatrise NIfTI load failed %s: %s', heatrise_files(1).name, ME.message);
                end
            else
                warning('Heatrise NIfTI not found: %s', heatrise_pat);
            end
        end
    end
end

%% Plot Isppa at post-skull voxel
plot_metric_with_patches(isppa_postskull, xlabels, 'Isppa at post-skull voxel (W/cm^2)', ...
    'sub-002_maxIsppa_skull_brain', [0 1], scatter_levels, scatter_colors, pn);

%% Plot riseT at post-skull voxel
plot_metric_with_patches(heatrise_postskull, xlabels, 'rise T brain (deg C)', ...
    'sub-002_riseT_skull_brain', [0 1.4], scatter_levels, scatter_colors, pn);

%% Skull-properties-only plots (analogous to benchmark2 skull-only section)
n_props       = length(property_labels);
skull_idx     = n_props+1 : 2*n_props;
xlabels_skull = cellfun(@(sym) [sym '_{s}'], property_symbols, 'UniformOutput', false);

plot_metric_with_patches(isppa_postskull(:, skull_idx), xlabels_skull, ...
    'Isppa at post-skull voxel (W/cm^2)', ...
    'sub-002_maxIsppa_skull_brain_skullonly', [0 1], scatter_levels, scatter_colors, pn);

plot_metric_with_patches(heatrise_postskull(:, skull_idx), xlabels_skull, ...
    'rise T brain (deg C)', ...
    'sub-002_riseT_skull_brain_skullonly', [0 1.4], scatter_levels, scatter_colors, pn);
