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
postskull_vox = [35, 35, 58]; % fixed post-skull voxel (x,y,z)

% Define tissues and suffixes
tissue_labels   = {'brain', 'skull'};
tissue_suffixes = {'_b', '_s'};

% Property definitions (must match simulation script)
property_labels   = {'sound_speed', 'density', 'alpha_coeff', ...
                     'thermal_conductivity', 'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};

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
        xlabels{i_cond} = [tissue_suf(2:end) '_' prop_label];

        for sw = sweep_values
            affix    = [tissue_suf prop_suf num2str(sw)];
            csv_file = fullfile(subdir_path, ['sub-002_phantom_output_table' affix '.csv']);

            if ~isfile(csv_file)
                warning('CSV not found: %s', csv_file);
                continue
            end

            T_csv = readtable(csv_file);

            % riseT at the post-skull location from the output table
            if ismember('riseT_brain', T_csv.Properties.VariableNames)
                heatrise_postskull(sw, i_cond) = T_csv.riseT_brain(1);
            else
                warning('riseT_brain missing in %s', csv_file);
            end

            % Isppa at the fixed post-skull voxel from the intensity NIfTI
            nifti_pat   = fullfile(subdir_path, sprintf('*final_intensity_orig_coord%s%s%d.nii.gz', tissue_suf, prop_suf, sw));
            nifti_files = dir(nifti_pat);
            if ~isempty(nifti_files)
                try
                    V = niftiread(fullfile(nifti_files(1).folder, nifti_files(1).name));
                    voxel_val = V(postskull_vox(1), postskull_vox(2), postskull_vox(3));
                    isppa_postskull(sw, i_cond) = voxel_val;
                    fprintf('Extracted %s sw%d (%s %s): voxel(%d,%d,%d)=%.4f\n', ...
                        subject_name, sw, tissue_labels{t}, prop_label, ...
                        postskull_vox(1), postskull_vox(2), postskull_vox(3), voxel_val);
                catch ME
                    warning('NIfTI load failed %s: %s', nifti_files(1).name, ME.message);
                end
            else
                warning('Intensity NIfTI not found for pattern: %s', nifti_pat);
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
