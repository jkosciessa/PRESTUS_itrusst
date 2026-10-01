% Extract focus-voxel metrics for sub-001 (brain-only benchmark).
% Reads Isppa at the empirical brain focus from the intensity NIfTI and
% riseT_brain from the output-table CSV (no heatrise NIfTI is produced).

% get root path (script must be run from file location)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'))
rootpath = pwd;

pn.figures = fullfile(rootpath, 'figures');
pn.data    = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a', 'sub-001');
pn.code    = fullfile(rootpath, 'code'); addpath(fullfile(pn.code, 'functions'));

% Define tissues and suffixes (sub-001 is brain-only benchmark)
tissue_labels   = {'brain'};
tissue_suffixes = {'_b'};

% Property definitions (must match simulation script)
property_labels   = {'sound_speed', 'density', 'alpha_coeff', ...
                     'thermal_conductivity', 'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};
property_symbols  = {'c', '\rho', '\alpha', '\kappa', 'c_p', '\omega', '\beta'};

sweep_values   = 1:2;
scatter_levels = [1, 2];
scatter_colors = {[1 1 1], [0 0 0]};

n_conds = length(tissue_labels) * length(property_labels);
isppa_focus   = nan(length(sweep_values), n_conds);
heatrise_focus = nan(length(sweep_values), n_conds);
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
            csv_file = fullfile(pn.data, ['sub-001_phantom' affix '.csv']);

            if ~isfile(csv_file)
                warning('CSV not found: %s', csv_file);
                continue
            end

            T_csv = readtable(csv_file);

            if ~ismember('Ix_brain_vox', T_csv.Properties.VariableNames)
                warning('Ix_brain_vox missing in %s', csv_file);
                continue
            end
            ix = T_csv.Ix_brain_vox(1);
            iz = T_csv.Iz_brain_vox(1);

            % Read Isppa at the empirical brain focus voxel from intensity NIfTI
            intensity_pat   = fullfile(pn.data, 'nii', sprintf('*phantom_T1w%s%s%d_intensity.nii.gz', tissue_suf, prop_suf, sw));
            intensity_files = dir(intensity_pat);
            if ~isempty(intensity_files)
                try
                    V = niftiread(fullfile(intensity_files(1).folder, intensity_files(1).name));
                    isppa_focus(i_sw, i_cond) = V(ix, iz);
                    fprintf('sub-001 sw%d (%s %s): Isppa voxel(%d,%d)=%.4f\n', ...
                        sw, tissue_labels{t}, prop_label, ix, iz, V(ix,iz));
                catch ME
                    warning('Intensity NIfTI load failed %s: %s', intensity_files(1).name, ME.message);
                end
            else
                warning('Intensity NIfTI not found: %s', intensity_pat);
            end

            % Read heatrise at the empirical brain focus voxel from heatrise NIfTI
            heatrise_pat   = fullfile(pn.data, 'nii', sprintf('*phantom_T1w%s%s%d_heatrise.nii.gz', tissue_suf, prop_suf, sw));
            heatrise_files = dir(heatrise_pat);
            if ~isempty(heatrise_files)
                try
                    Vh = niftiread(fullfile(heatrise_files(1).folder, heatrise_files(1).name));
                    heatrise_focus(i_sw, i_cond) = Vh(ix, iz);
                    fprintf('sub-001 sw%d (%s %s): heatrise voxel(%d,%d)=%.4f\n', ...
                        sw, tissue_labels{t}, prop_label, ix, iz, Vh(ix,iz));
                catch ME
                    warning('Heatrise NIfTI load failed %s: %s', heatrise_files(1).name, ME.message);
                end
            else
                warning('Heatrise NIfTI not found: %s', heatrise_pat);
            end
        end
    end
end

%% Plot Isppa at empirical focus
plot_metric_with_patches(isppa_focus, xlabels, 'Isppa at focus (W/cm^2)', ...
    'sub-001_maxIsppa_focus', [0 60], scatter_levels, scatter_colors, pn);

%% Plot riseT at empirical focus
plot_metric_with_patches(heatrise_focus, xlabels, 'rise T brain (deg C)', ...
    'sub-001_riseT_focus', [0 1.5], scatter_levels, scatter_colors, pn);
