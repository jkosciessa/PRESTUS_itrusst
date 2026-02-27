clear all; close all; clc;

user = 'julkos'

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code'));

% use most recent PRESTUS update (merge)
pn.tuSIM = fullfile(rootpath, '..', 'merge', 'tools', 'PRESTUS'); addpath(pn.tuSIM);
pn.tuSIM_fun = fullfile(pn.tuSIM, 'functions'); addpath(genpath(pn.tuSIM_fun));
pn.tuSIM_profile = fullfile(pn.tuSIM, 'acoustic_profiling'); addpath(genpath(pn.tuSIM_profile));
pn.tuSIM_tools = fullfile(pn.tuSIM, 'toolboxes'); addpath(genpath(pn.tuSIM_tools));
pn.kwave = fullfile(pn.tuSIM_tools, 'k-wave', 'k-Wave'); addpath(pn.kwave);
if strcmp(user, 'julkos')
    pn.simnibs = fullfile('/home', 'neuromod', 'julkos', '.conda', 'envs', 'simnibs_env', 'bin');
else
    warning('User not specified...may be missing simnibs path');
end
pn.minimize = fullfile(pn.tuSIM_tools, 'FEX-minimize'); addpath(pn.minimize);
pn.configs = fullfile(rootpath, 'data', 'configs');
pn.data_path = fullfile(rootpath, 'data', 'bids');
pn.data_seg = fullfile(rootpath, 'data', 'simnibs');
pn.data = fullfile(rootpath, 'data');
pn.nifti = (fullfile(rootpath, 'tools', 'nifti_toolbox')); addpath(pn.nifti);

%% define variables to loop

setup_list = {['itrusst_protocol1']};

all_subjects = 003;
intensities = 60;
copy_positions = 1;

for subject_id = all_subjects

    for i_setup = 1:length(setup_list)
    for i_intensity = 1:length(intensities)
        transducer_name = setup_list{i_setup};
        desired_intensity = intensities(i_intensity);

        %% load parameter file and adjust paths if necessary
        
        cd(fullfile(rootpath, 'data', 'configs'))
        parameters = load_parameters(['config_',transducer_name,'.yaml']);

        if ismac
            % reset paths to local mac install; NOTE: horrible not to have this be dynamic
            parameters.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
        else
            parameters.simnibs_bin_path = pn.simnibs;
        end
        pn.data_sims = fullfile(rootpath, 'data', 'tussim', [transducer_name]);
            if ~exist(pn.data_sims); mkdir(pn.data_sims); end
        parameters.ld_library_path ="/opt/gcc/7.2.0/lib64";
        parameters.data_path = pn.data_seg; % use simnibs folder
        parameters.seg_path = pn.data_seg;
        parameters.t1_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
        parameters.t2_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
        parameters.sim_path = pn.data_sims;
        parameters.paths_to_add = {pn.kwave, pn.minimize};
        pn.outputs_folder = fullfile(parameters.sim_path,sprintf('sub-%03d', subject_id));
        if ~exist(pn.outputs_folder); mkdir(pn.outputs_folder); end
        parameters.usepseudoCT = 0;

        calibrate = 0;
        if calibrate == 1
            [opt_source_amp, opt_source_phase_deg, opt_source_phase_rad] = ...
                itrusst_calibrate(pn, parameters, subject_id, desired_intensity, transducer_name);
        end

        %% run soft-tissue simulation

        cd(fullfile(rootpath, 'data', 'configs'))
        parameters = load_parameters(['config_',transducer_name,'.yaml']);

        parameters.ld_library_path ="/opt/gcc/7.2.0/lib64";
        parameters.data_path = pn.data_seg;
        parameters.seg_path = pn.data_seg;
        parameters.t1_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
        parameters.t2_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
        parameters.sim_path = pn.data_sims;
        parameters.paths_to_add = {pn.kwave, pn.minimize};

        parameters.simulation_medium = 'phantom'; % use default grid in combo with (final) tissue mask
        parameters.interactive = 0;
        parameters.overwrite_files = 'always';
        parameters.overwrite_simnibs = 0;
        parameters.heatingvideo = 0;

        if calibrate == 1
            parameters.transducer.source_amp = opt_source_amp; % use calibrated input intensity
            parameters.transducer.source_phase_deg = opt_source_phase_deg;
            parameters.transducer.source_phase_rad = opt_source_phase_rad;
        end

        parameters.transducer.trans_pos = [3, 35];
        parameters.transducer.focus_pos = [64, 35];

        parameters.run_source_setup = 1;
        parameters.run_acoustic_sims = 1;
        parameters.run_heating_sims = 1;
        parameters.run_posthoc_water_sims = 1;
        parameters.axisymmetric = 1;

        %% iterate over tissue parameters
        % benchmark sub-001: brain
        % benchmark sub-002: brain & skull
        % [x] benchmark sub-003: brain & skull trabecular & skull cortical 

        % iterate over multiple parameters:
        % [1] the tissue labels "skull-trabecular" and "skull-cortical", create either a suffix "_st_" or "_sc_"
        % [2] for each of these tissues, there are multiple property labels: "sound_speed", "density", "alpha_0_true", "alpha_power_true", "thermal_conductivity", "specific_heat_capacity", create a short suffix 
        % 
        % For each of these combinations:
        % - specify a "parameters.results_filename_affix" as a combination of these suffixes
        % - vary upper and lower parameters based on tissue_parameters_multiskull

        % Define tissue labels and suffixes
        tissue_labels = {'skull_trabecular', 'skull_cortical'};
        tissue_suffixes = {'_st', '_sc'};
        
        % Define property labels and short suffixes
        property_labels = {'sound_speed', 'density', 'alpha_0_true','thermal_conductivity', ...
                           'specific_heat_capacity', 'perfusion', ...
                           'absorption_fraction'};
        property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};

        % define loops
        magnitude = {'bottom', 'top'};
        
        % load ranges
        opts = detectImportOptions(fullfile(pn.data, 'tissue_parameters_multiskull.csv'));
        params_table = readtable(fullfile(pn.data, 'tissue_parameters_multiskull.csv'), opts); clear opts;

        % Loop over tissues
        for t = 1:length(tissue_labels)
        
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
                            parameters.medium.(tissue_labels{tissue_bl}).(property_labels{pp}) = ...
                                params_table.(property_labels{pp})(find(strcmp(params_table.tissue,tissue_labels{tissue_bl}) & strcmp(params_table.level,'mid')));
                        end
                    end
                    
                    % Vary only the current property to BOTTOM/TOP
                    parameters.medium.(tissue).(property) = ...
                        params_table.(property_labels{p})(find(strcmp(params_table.tissue,tissue) & strcmp(params_table.level,magnitude{v})));

                    % Create filename affix
                    parameters.results_filename_affix = ...
                        [tissue_suffix prop_suffix num2str(v)];
        
                    % Run the simulation
                    parameters.code_type = 'matlab_gpu';
                    parameters.hpc_submit_medium = 'matlab';

                    switch parameters.hpc_submit_medium
                        case 'slurm'
                            single_subject_pipeline_with_slurm(subject_id, parameters, 0, '01:00:00', 20);
                        case 'matlab'
                            single_subject_pipeline(subject_id, parameters);
                    end

                end
            end
        end

        % end loop

    end
    end
end
