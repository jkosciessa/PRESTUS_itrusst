clear all; close all; clc;

user = 'julkos' % julkos/marwim/sebrei

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code'));

pn.prestus = fullfile(rootpath, '..', 'v05', 'tools', 'PRESTUS'); addpath(genpath(pn.prestus));
pn.prestus_fun = fullfile(pn.prestus, 'functions'); addpath(genpath(pn.prestus_fun));
pn.prestus_profile = fullfile(pn.prestus, 'acoustic_profiling'); addpath(genpath(pn.prestus_profile));
pn.prestus_tools = fullfile(pn.prestus, 'toolboxes'); addpath(genpath(pn.prestus_tools));
pn.kwave = fullfile(pn.prestus_tools, 'k-wave', 'k-Wave'); addpath(pn.kwave);
if strcmp(user, 'julkos')
    pn.simnibs = fullfile('/home', 'neuromod', 'julkos', '.conda', 'envs', 'simnibs_env', 'bin');
elseif strcmp(user, 'marwim')
    pn.simnibs = fullfile('/home', 'neuromod', 'marwim', 'miniconda3', 'envs', 'simnibs_env', 'bin');
elseif strcmp(user, 'sebrei')
    pn.simnibs = fullfile('/home', 'neuromod', 'sebrei', '.conda', 'envs', 'simnibs_env', 'bin');
end
pn.minimize = fullfile(pn.prestus_tools, 'FEX-minimize'); addpath(pn.minimize);
pn.configs = fullfile(rootpath, 'data', 'configs');
pn.data_path = fullfile(rootpath, 'data', 'bids');
pn.path_seg = fullfile(rootpath, 'data', 'simnibs');
pn.data = fullfile(rootpath, 'data');
pn.nifti = (fullfile(rootpath, 'tools', 'nifti_toolbox')); addpath(pn.nifti);

%% define variables to loop

setup_list = {['itrusst_protocol1']};

all_subjects = 301:312;
intensities = 60;
copy_positions = 1;
calibrate = 0;

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
            parameters.startup.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
        else
            parameters.startup.simnibs_bin_path = pn.simnibs;
        end
        pn.path_sim = fullfile(rootpath, 'data', 'tussim', [transducer_name]);
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

        if calibrate == 1
            [opt_source_amp, opt_source_phase_deg, opt_source_phase_rad] = ...
                itrusst_calibrate(pn, parameters, subject_id, desired_intensity, transducer_name);
        end

        %% run soft-tissue simulation

        cd(fullfile(rootpath, 'data', 'configs'))
        parameters = load_parameters(['config_',transducer_name,'.yaml']);

        parameters.hpc.ld_library_path ="/opt/gcc/7.2.0/lib64";
        parameters.path.anat = pn.path_seg;
        parameters.path.seg = pn.path_seg;
        parameters.path.t1_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
        parameters.path.t2_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
        parameters.path.sim = pn.path_sim;
        parameters.startup.paths_to_add = {pn.kwave, pn.minimize};

        parameters.simulation.medium = 'phantom'; % use default grid in combo with (final) tissue mask
        parameters.simulation.interactive = 0;
        parameters.io.overwrite_files = 'always';
        parameters.io.overwrite_simnibs = 0;
        parameters.io.save_heatingvideo = 0;

        if calibrate == 1
            parameters.transducer.annular.elem_amp = opt_source_amp;
            parameters.transducer.annular.elem_phase_deg = opt_source_phase_deg;
            parameters.transducer.annular.elem_phase_rad = opt_source_phase_rad;
        end

        %% iterate over tissue parameters

        % Define tissue labels and suffixes
        tissue_labels = {'brain', 'skull'};
        tissue_suffixes = {'_b', '_s'};
        
        % Define property labels and short suffixes
        property_labels = {'sound_speed', 'density', 'alpha_coeff','thermal_conductivity', ...
                           'specific_heat_capacity', 'perfusion', ...
                           'absorption_fraction'};
        property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};

        % define loops
        magnitude = {'mid'};
        
        % load ranges
        opts = detectImportOptions(fullfile(pn.data, 'tissue_parameters.csv'));
        params_table = readtable(fullfile(pn.data, 'tissue_parameters.csv'), opts); clear opts;

        % Loop over tissues
        for t = 1:length(tissue_labels)
            tissue = tissue_labels{t};
            tissue_suffix = tissue_suffixes{t};

            % Set all properties to 100% (baseline)
            for pp = 1:length(property_labels)
                parameters.medium_properties.(tissue).(property_labels{pp}) = ...
                    params_table.(property_labels{pp})(find(strcmp(params_table.tissue,tissue) & strcmp(params_table.level,'mid')));
            end
        end 
        
        % Create filename affix
        parameters.io.output_affix = "";

        % Run the simulation
        parameters.simulation.code_type = 'cpp_cpu';
        parameters.simulation.precision = 'single';

        parameters.grid.resolution_mm = 0.5;
        parameters.grid.default_dims = [70, 120]/parameters.grid.resolution_mm;
        parameters.grid.pml_size = [10]/parameters.grid.resolution_mm;

        parameters.modules.run_heating_sims = 0;

        parameters.subject_id = subject_id;
        prestus_pipeline_start(parameters);

        % end loop

    end
    end
end
