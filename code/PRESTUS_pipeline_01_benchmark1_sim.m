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
pn = setup_itrusst_paths(rootpath, user);

%% define variables to loop

setup_list = {
    'medium_conservative'; ...
    'medium_liberal'; ...
    'medium_mid'; ...
    'medium_BabelBrain'; ...
    'medium_kPlan'; ...
    'medium_PRESTUS'};

all_subjects = [001:002]; %602
intensities = 60;
copy_positions = 1;
calibrate = 1;

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
            parameters.startup.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
        else
            parameters.startup.simnibs_bin_path = pn.simnibs;
        end
        pn.path_sim = fullfile(rootpath, 'data', 'tussim', [setup_name]);
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
                itrusst_calibrate(pn, parameters, subject_id, desired_intensity, setup_name);
        end

        %% run simulation

        cd(fullfile(rootpath, 'data', 'configs'));
        
        parameters = load_parameters(base_config, pwd, ['config_',setup_name,'.yaml'], pwd);

        parameters.hpc.ld_library_path ="/opt/gcc/7.2.0/lib64";
        parameters.path.anat = pn.path_seg;
        parameters.path.seg = pn.path_seg;
        parameters.path.t1_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
        parameters.path.t2_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
        parameters.path.sim = pn.path_sim;
        parameters.startup.paths_to_add = {pn.kwave, pn.minimize};

        parameters.simulation.interactive = 0;
        parameters.io.overwrite_files = 'always';
        parameters.io.overwrite_simnibs = 0;
        parameters.io.save_heatingvideo = 0;

        if calibrate == 1
            parameters.transducer.annular.elem_amp = opt_source_amp; % use calibrated input intensity
            parameters.transducer.annular.elem_phase_deg = opt_source_phase_deg;
            parameters.transducer.annular.elem_phase_rad = opt_source_phase_rad;
        end
        
        if subject_id == 602 % for 3D run layered with C++ GPU code
            parameters.simulation.code_type = 'cpp_gpu';
            parameters.simulation.medium = 'layered';
            parameters.grid.resolution_mm = 0.75;
            parameters.grid.axisymmetric = 0;
            jobtime = '02:00:00';
            jobram = 20;
        else
            parameters.simulation.code_type = 'matlab_cpu';
            parameters.simulation.medium = 'phantom';
            parameters.grid.axisymmetric = 1;
            jobtime = '00:30:00';
            jobram = 8;
        end

        % specify GPU setup if requested
        if contains(parameters.simulation.code_type, 'gpu')
            parameters.hpc.gpu = "gpu:1";
            parameters.hpc.partition = "gpu40g";
        end

        % Run the simulation
        parameters.subject_id = subject_id;
        prestus_pipeline_start(parameters);

    end
    end
end
