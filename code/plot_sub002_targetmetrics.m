% Extract focus-voxel metrics for sub-002 (brain + skull benchmark).
% Reads Isppa at the empirical brain focus from the intensity NIfTI and
% riseT_brain from the output-table CSV (no heatrise NIfTI is produced).

% get root path (script must be run from file location)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.figures = fullfile(rootpath, 'figures');
pn.data    = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a', 'sub-002');
pn.code    = fullfile(rootpath, 'code'); addpath(pn.code);

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
isppa_focus    = nan(length(sweep_values), n_conds);
heatrise_focus = nan(length(sweep_values), n_conds);
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
            csv_file = fullfile(pn.data, ['sub-002_phantom_output_table' affix '.csv']);

            if ~isfile(csv_file)
                warning('CSV not found: %s', csv_file);
                continue
            end

            T_csv = readtable(csv_file);

            % Read riseT_brain directly from the output table
            if ismember('riseT_brain', T_csv.Properties.VariableNames)
                heatrise_focus(sw, i_cond) = T_csv.riseT_brain(1);
            else
                warning('riseT_brain missing in %s', csv_file);
            end

            % Read Isppa at the empirical brain focus voxel from intensity NIfTI
            nifti_pat   = fullfile(pn.data, sprintf('*final_intensity_orig_coord%s%s%d.nii.gz', tissue_suf, prop_suf, sw));
            nifti_files = dir(nifti_pat);
            if ~isempty(nifti_files) && ...
               ismember('Ix_brain_vox', T_csv.Properties.VariableNames)
                try
                    V = niftiread(fullfile(nifti_files(1).folder, nifti_files(1).name));
                    ix = T_csv.Ix_brain_vox(1);
                    iy = T_csv.Iy_brain_vox(1);
                    iz = T_csv.Iz_brain_vox(1);
                    isppa_focus(sw, i_cond) = V(ix, iy, iz);
                    fprintf('sub-002 sw%d (%s %s): voxel(%d,%d,%d)=%.4f\n', ...
                        sw, tissue_labels{t}, prop_label, ix, iy, iz, V(ix,iy,iz));
                catch ME
                    warning('NIfTI load failed %s: %s', nifti_files(1).name, ME.message);
                end
            elseif isempty(nifti_files)
                warning('Intensity NIfTI not found for pattern: %s', nifti_pat);
            end
        end
    end
end

%% Plot Isppa at empirical focus
plot_metric_with_patches(isppa_focus, xlabels, 'Isppa at focus (W/cm^2)', ...
    'sub-002_maxIsppa_focus', [0 10], scatter_levels, scatter_colors, pn);

%% Plot riseT at empirical focus
plot_metric_with_patches(heatrise_focus, xlabels, 'rise T brain (deg C)', ...
    'sub-002_riseT_focus', [0 0.25], scatter_levels, scatter_colors, pn);
