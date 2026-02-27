% This script executes PRESTUS stimulations for multiple config files
% differing in their parameter assumptions.

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
pn.nifti = (fullfile(rootpath, 'tools', 'nifti_toolbox')); addpath(pn.nifti);

%% define variables to loop

setup_list = {
    'medium_conservative'; ...
    'medium_liberal'; ...
    'medium_mid'; ...
    'medium_BabelBrain'; ...
    'medium_kPlan'; ...
    'medium_PRESTUS'};

all_subjects = 602;% [001:002];
intensities = 60;
copy_positions = 1;
calibrate = 0;

for subject_id = all_subjects

    for i_setup = 1:length(setup_list)
    for i_intensity = 1:length(intensities)
        setup_name = setup_list{i_setup};
        desired_intensity = intensities(i_intensity);

        %% load parameter file and adjust paths if necessary
        
        cd(fullfile(rootpath, 'data', 'configs'));

        if subject_id == 602
            base_config = 'config_itrusst_protocol1_602.yaml';
        else
            base_config = 'config_itrusst_protocol1.yaml';
        end

        parameters = load_parameters(base_config, pwd, ['config_',setup_name,'.yaml'], pwd);

        if ismac
            % reset paths to local mac install; NOTE: horrible not to have this be dynamic
            parameters.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
        else
            parameters.simnibs_bin_path = pn.simnibs;
        end
        pn.data_sims = fullfile(rootpath, 'data', 'tussim', [setup_name]);
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

        if calibrate == 1
            [opt_source_amp, opt_source_phase_deg, opt_source_phase_rad] = ...
                itrusst_calibrate(pn, parameters, subject_id, desired_intensity, setup_name);
        end

        %% run simulation

        parameters = load_parameters(base_config, pwd, ['config_',setup_name,'.yaml'], pwd);

        parameters.ld_library_path ="/opt/gcc/7.2.0/lib64";
        parameters.data_path = pn.data_seg;
        parameters.seg_path = pn.data_seg;
        parameters.t1_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
        parameters.t2_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
        parameters.sim_path = pn.data_sims;
        parameters.paths_to_add = {pn.kwave, pn.minimize};

        parameters.interactive = 0;
        parameters.overwrite_files = 'always';
        parameters.overwrite_simnibs = 0;
        parameters.heatingvideo = 0;

        if calibrate == 1
            parameters.transducer.source_amp = opt_source_amp; % use calibrated input intensity
            parameters.transducer.source_phase_deg = opt_source_phase_deg;
            parameters.transducer.source_phase_rad = opt_source_phase_rad;
        end
        
        if subject_id == 602 % for 3D run layered with C++ GPU code
            parameters.code_type = 'cpp_gpu';
            parameters.simulation_medium = 'layered';
            parameters.grid_step_mm = 0.75;
            parameters.axisymmetry = 0;
            jobtime = '02:00:00';
            jobram = 20;
        else
            parameters.code_type = 'matlab_cpu';
            parameters.simulation_medium = 'phantom';
            parameters.axisymmetry = 1;
            jobtime = '00:30:00';
            jobram = 8;
        end

        % specify GPU setup if requested
        if contains(parameters.code_type, 'gpu')
            parameters.hpc_gpu = "gpu:1";
            parameters.hpc_partition = "gpu40g";
        end

        parameters.hpc_submit_medium = 'slurm';
        
        % Run the simulation
        switch parameters.hpc_submit_medium
            case 'slurm'
                single_subject_pipeline_with_slurm(subject_id, parameters, 0, jobtime, jobram);
            case 'matlab'
                single_subject_pipeline(subject_id, parameters);
        end

    end
    end
end
