clear all; close all; clc;

user = 'julkos'

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
rootpath = fileparts(pathstr);

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, user);

%% define variables to loop

% 500 kHz @ 60 W/cm2 (baseline)
variant_names     = {'itrusst_protocol1a'};
variant_configs   = {'itrusst_protocol'};
variant_intensity = [60];

all_subjects = [001, 002, 303];

for subject_id = all_subjects

    for i_variant = 1:length(variant_names)
        setup_name        = variant_configs{i_variant};
        variant_name      = variant_names{i_variant};
        desired_intensity = variant_intensity(i_variant);

        %% load parameter file and adjust paths if necessary
        
        cd(fullfile(rootpath, 'data', 'configs'))
        parameters = load_parameters(['config_',setup_name,'.yaml']);

        if ismac
            % reset paths to local mac install; NOTE: horrible not to have this be dynamic
            parameters.startup.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
        else
            parameters.startup.simnibs_bin_path = pn.simnibs;
        end
        pn.path_sim = fullfile(rootpath, 'data', 'tussim', variant_name);
            if ~exist(pn.path_sim); mkdir(pn.path_sim); end
        parameters.hpc.ld_library_path ="/opt/gcc/7.2.0/lib64";
        parameters.path.anat = pn.path_seg; % use simnibs folder
        parameters.path.seg = pn.path_seg;
        parameters.path.t1_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
        parameters.path.t2_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
        parameters.path.sim = pn.path_sim;
        parameters.startup.paths_to_add = {pn.kwave, pn.minimize};
        pn.outputs_folder = fullfile(pn.path_sim,sprintf('sub-%03d', subject_id));
        if ~exist(pn.outputs_folder); mkdir(pn.outputs_folder); end
        parameters.pct.enabled = 0;

        %% run layered phantom simulation

        cd(fullfile(rootpath, 'data','configs'));
        parameters = load_parameters(['config_',setup_name,'.yaml']);

        parameters.hpc.ld_library_path ="/opt/gcc/7.2.0/lib64";
        parameters.path.anat = pn.path_seg;
        parameters.path.seg = pn.path_seg;
        parameters.path.t1_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
        parameters.path.t2_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
        parameters.path.sim = pn.path_sim;
        parameters.startup.paths_to_add = {pn.kwave, pn.minimize};

        parameters.simulation.medium = 'phantom'; % use default grid in combo with (final) tissue mask
        parameters.simulation.interactive = 0;
        parameters.io.overwrite_files = 'never';
        parameters.io.overwrite_simnibs = 0;
        parameters.io.save_heatingvideo = 0;
        parameters.io.save_matrices = 0;

        %% iterate over tissue parameters
        % benchmark sub-001: brain
        % benchmark sub-002: brain & skull

        % iterate over multiple parameters:
        % [1] the tissue labels "brain" and "skull", create either a suffix "_b_" or "_s_"
        % [2] for each of these tissues, there are multiple property labels: "sound_speed", "density", "alpha_0_true", "alpha_power_true", "thermal_conductivity", "specific_heat_capacity", create a short suffix 
        % 
        % For each of these combinations:
        % - specify a "parameters.results_filename_affix" as a combination of these suffixes
        % - vary "parameters_baseline.medium.<<brain>>.<<sound_speed>>" from 10% to 200% in steps of 10 %, set the other parameters to their respective 100% value.

        % Define tissue labels and suffixes
        tissue_labels = {'brain', 'skull'};
        tissue_suffixes = {'_b', '_s'};
        
        % Define property labels and short suffixes
        property_labels = {'sound_speed', 'density', 'alpha_coeff','thermal_conductivity', ...
                           'specific_heat_capacity', 'perfusion', ...
                           'absorption_fraction'};
        property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};

        % define loops
        magnitude = {'bottom', 'top'};
        
        % load ranges
        opts = detectImportOptions(fullfile(pn.data, 'tissue_parameters.csv'));
        params_table = readtable(fullfile(pn.data, 'tissue_parameters.csv'), opts); clear opts;

        % Loop over tissues
        for t = 1:length(tissue_labels)
        
            % no skull exists for benchmark 1
            if subject_id==1 && t==2
                continue;
            end
        
            tissue = tissue_labels{t};
            tissue_suffix = tissue_suffixes{t};
            
            % Loop over properties
            for p = 1:length(property_labels)
                property = property_labels{p};
                prop_suffix = property_suffixes{p};
                                
                % Loop over sweep values for the current property
                for v = 1:numel(magnitude)
                    % Set all properties to mid setting (baseline)
                    for tissue_bl = 1:length(tissue_labels)
                        for pp = 1:length(property_labels)
                            parameters.medium_properties.(tissue_labels{tissue_bl}).(property_labels{pp}) = ...
                                params_table.(property_labels{pp})(find(strcmp(params_table.tissue,tissue_labels{tissue_bl}) & strcmp(params_table.level,'mid')));
                        end
                    end
                    
                    % Vary only the current property to BOTTOM/TOP
                    parameters.medium_properties.(tissue).(property) = ...
                        params_table.(property_labels{p})(find(strcmp(params_table.tissue,tissue) & strcmp(params_table.level,magnitude{v})));

                    % Create filename affix
                    parameters.io.output_affix = ...
                        [tissue_suffix prop_suffix num2str(v)];

                    % Run the simulation
                    parameters.simulation.code_type = 'cpp_cpu'; % use CPU here to avoid parallel GPU calls
                    parameters.simulation.precision = 'single';

                    parameters.grid.resolution_mm = 0.5;
                    parameters.grid.default_dims = [70, 120]/parameters.grid.resolution_mm;
                    parameters.grid.pml_size = [20]/parameters.grid.resolution_mm; % 20 mm PML, matches primary pipeline configs
                    parameters.hpc.timelimit = '04:00:00';
                    parameters.grid.axisymmetric = 1;
  
                    parameters.subject_id = subject_id;
                    prestus_pipeline_start(parameters);
                end
            end
        end

        % end loop

    end
end
