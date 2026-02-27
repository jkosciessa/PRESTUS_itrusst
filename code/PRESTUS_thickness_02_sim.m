clear all; close all; clc;

user = 'julkos' % julkos/marwim/sebrei

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code'));

pn.tuSIM = fullfile(rootpath, 'tools', 'PRESTUS'); addpath(pn.tuSIM);
pn.tuSIM_fun = fullfile(pn.tuSIM, 'functions'); addpath(pn.tuSIM_fun);
pn.tuSIM_tools = fullfile(pn.tuSIM, 'toolboxes'); addpath(genpath(pn.tuSIM_tools));
pn.kwave = fullfile(pn.tuSIM_tools, 'k-wave', 'k-Wave'); addpath(pn.kwave);
if strcmp(user, 'julkos')
    pn.simnibs = fullfile('/home', 'neuromod', 'julkos', '.conda', 'envs', 'simnibs_env', 'bin');
elseif strcmp(user, 'marwim')
    pn.simnibs = fullfile('/home', 'neuromod', 'marwim', 'miniconda3', 'envs', 'simnibs_env', 'bin');
elseif strcmp(user, 'sebrei')
    pn.simnibs = fullfile('/home', 'neuromod', 'sebrei', '.conda', 'envs', 'simnibs_env', 'bin');
end
pn.minimize = fullfile(pn.tuSIM_tools, 'FEX-minimize'); addpath(pn.minimize);
pn.configs = fullfile(rootpath, 'data', 'configs');
pn.data_path = fullfile(rootpath, 'data', 'bids');
pn.data_seg = fullfile(rootpath, 'data', 'simnibs');
pn.data = fullfile(rootpath, 'data');
pn.nifti = (fullfile(rootpath, 'tools', 'nifti_toolbox')); addpath(pn.nifti);

%% define variables to loop

setup_list = {['itrusst_protocol1']};

all_subjects = 301:312;
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

        %% run soft-tissue simulation

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

        %% iterate over tissue parameters

        % Define tissue labels and suffixes
        tissue_labels = {'brain', 'skull'};
        tissue_suffixes = {'_b', '_s'};
        
        % Define property labels and short suffixes
        property_labels = {'sound_speed', 'density', 'alpha_0_true','thermal_conductivity', ...
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
                parameters.medium.(tissue).(property_labels{pp}) = ...
                    params_table.(property_labels{pp})(find(strcmp(params_table.tissue,tissue) & strcmp(params_table.level,'mid')));
            end
        end 
        
        % Create filename affix
        parameters.results_filename_affix = "";

        % Run the simulation
        parameters.code_type = 'matlab_cpu';
        single_subject_pipeline_with_slurm(subject_id, parameters, 0, '01:00:00', 16);

        % end loop

    end
    end
end
