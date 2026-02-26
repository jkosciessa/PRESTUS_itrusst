%% Fixed plotting script - CORRECTED table init + FULL nested loops
% Uses string array for VariableTypes (MATLAB R2019b+ compatible)

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.figures = fullfile(rootpath, 'figures');
pn.data = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1');
pn.code = fullfile(rootpath, 'code'); addpath(pn.code);

% Define tissues and suffixes
tissue_labels = {'brain', 'skull'};
tissue_suffixes = {'_b', '_s'};

% Define properties and suffixes
property_labels = {'soundSpeed', 'density', 'attenuation', ...
                   'thermConductivity', 'heatCapacity', 'perfusion', 'absorption'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};

% Define sweep values (2 levels)
sweep_values = 1:2;

% FIXED: Correct table initialization (string array for types)
var_names = {'subdir', 'sweep_level', 'tissue', 'property', 'property_value', ...
             'riseT_brain', 'max_Isppa_brain', 'isppa_postskull', 'heatrise_postskull'};
var_types = ["string", "double", "string", "string", "double", "double", "double", "double", "double"];
T = table('Size',[0 length(var_names)], 'VariableTypes', var_types, 'VariableNames', var_names);

v = 0;
subdir_path = fullfile(pn.data, 'sub-002');
    
for t = 1:length(tissue_labels)
    tissue_suf = tissue_suffixes{t};
    
    for p = 1:length(property_labels)
        v = v + 1;
        prop_label = property_labels{p};
        prop_suf = property_suffixes{p};
        
        for sw = sweep_values
            
            % Base categorical
            T.subdir(v) = subdirs(s_idx).name;
            T.sweep_level(v) = sw;
            T.tissue(v) = tissue_labels{t};
            T.property(v) = prop_label;
            
            % Property from CSV
            csv_pattern = fullfile(subdir_path, sprintf('*%s%s*%d*.csv', tissue_suf, prop_suf, sw));
            csv_files = dir(csv_pattern);
            prop_val = NaN;
            if ~isempty(csv_files)
                csv_file = fullfile(csv_files(1).folder, csv_files(1).name);
                try
                    T_csv = readtable(csv_file);
                    num_cols = varfun(@isnumeric, T_csv, 'OutputFormat', 'uniform');
                    if any(num_cols)
                        prop_val = T_csv{1, find(num_cols,1)};  % First numeric
                    end
                catch
                    warning('CSV read failed: %s', csv_file);
                end
            end
            T.property_value(v) = prop_val;
            
            % Priors (brain-specific, adapt paths/cols)
            T.riseT_brain(v) = NaN;
            T.max_Isppa_brain(v) = NaN;
            brain_csv_pat = fullfile(subdir_path, sprintf('*%s*_overview*%d.csv', tissue_suf, sw));
            brain_csvs = dir(brain_csv_pat);
            if ~isempty(brain_csvs)
                try
                    T_csv = readtable(fullfile(brain_csvs(1).folder, brain_csvs(1).name));
                    if ismember('riseT_brain', T_csv.Properties.VariableNames)
                        T.riseT_brain(v) = T_csv.riseT_brain(1);
                    end
                    if ismember('max_Isppa_brain', T_csv.Properties.VariableNames)
                        T.max_Isppa_brain(v) = T_csv.max_Isppa_brain(1);
                    end
                catch, end
            end
            
            isppa_val = NaN; heatrise_val = NaN;

            % NEW: Post-skull voxel (35,58) from BRAIN isppa NIfTI (sw-specific)
            nifti_pat = fullfile(subdir_path, sprintf('*final_isppa_orig_coord%s%s%d.nii.gz',tissue_suf, prop_suf, sw));
            nifti_files = dir(nifti_pat);
            if ~isempty(nifti_files)
                nifti_file = fullfile(nifti_files(1).folder, nifti_files(1).name);
                try
                    V = niftiread(nifti_file);
                    voxel_val = V(35,35,58);
                    isppa_val = voxel_val;
                    fprintf('Extracted %s sw%d (%s %s): voxel(35,58,%d)=%.4f\n', ...
                            subdirs(s_idx).name, sw, tissue_labels{t}, prop_label, midz, voxel_val);
                catch ME
                    warning('NIfTI load failed %s: %s', nifti_file, ME.message);
                end
            end

            % NEW: Post-skull voxel (35,58) from BRAIN heatrise NIfTI (sw-specific)
            nifti_pat = fullfile(subdir_path, sprintf('*final_heatrise_orig_coord%s%s%d.nii.gz',tissue_suf, prop_suf, sw));
            nifti_files = dir(nifti_pat);
            if ~isempty(nifti_files)
                nifti_file = fullfile(nifti_files(1).folder, nifti_files(1).name);
                try
                    V = niftiread(nifti_file);
                    voxel_val = V(35,35,58);
                    heatrise_val = voxel_val;
                    fprintf('Extracted %s sw%d (%s %s): voxel(35,58,%d)=%.4f\n', ...
                            subdirs(s_idx).name, sw, tissue_labels{t}, prop_label, midz, voxel_val);
                catch ME
                    warning('NIfTI load failed %s: %s', nifti_file, ME.message);
                end
            end

            isppa_postskull(sw, v) = isppa_val;
            heatrise_postskull(sw, v) = heatrise_val;
        end
    end
end

%% Plot max_Isppa_skull
plot_metric_with_patches(isppa_postskull, xlabels, 'max Isppa skull_brain', 'sub-002_maxIsppa_skull_brain', [0 1], scatter_levels, scatter_colors, pn);

%% Plot riseT_skull
plot_metric_with_patches(heatrise_postskull, xlabels, 'rise T skull_brain (deg C)', 'sub-002_riseT_skull_brain', [0 1.4], scatter_levels, scatter_colors, pn);
