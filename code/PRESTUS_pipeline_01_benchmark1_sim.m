% This script executes PRESTUS stimulations for multiple config files
% differing in their parameter assumptions.

clear all; close all; clc;

user = 'julkos' % julkos/marwim/sebrei

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, user);

%% define variables to loop

% literature setups: config_medium_heat_max/heat_min/trans_max/trans_min.yaml
% are regenerated here from data/tissue_parameters.csv + tissue_setups.csv
% (the reference medium is not part of this comparison)
itrusst_write_medium_configs(pn);

setup_list = {
    'medium_heat_max'; ...
    'medium_heat_min'; ...
    'medium_trans_max'; ...
    'medium_trans_min'; ...
    'medium_BabelBrain'; ...
    'medium_kPlan'; ...
    'medium_PRESTUS'};
% subset to submit, e.g. {'medium_kPlan'}; {} = all
run_setups = {};
if ~isempty(run_setups); setup_list = run_setups; end

code_type = 'cpp_gpu'; % gpu40g: ~3 min per job (A100, 2026-10-16); cpp_cpu: 70-105+ min

all_subjects = [101, 102];
intensities = 60;
copy_positions = 1;

for subject_id = all_subjects

    for i_setup = 1:length(setup_list)
    for i_intensity = 1:length(intensities)
        setup_name = setup_list{i_setup};
        desired_intensity = intensities(i_intensity);

        %% load parameter file and adjust paths if necessary
        
        cd(fullfile(rootpath, 'data', 'configs'));

        base_config = 'config_itrusst_protocol.yaml';

        parameters = load_parameters(base_config, pwd, ['config_',setup_name,'.yaml'], pwd);

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

        %% run simulation

        cd(fullfile(rootpath, 'data', 'configs'));
        
        parameters = load_parameters(base_config, pwd, ['config_',setup_name,'.yaml'], pwd);

        parameters.transducer.target_isppa_wcm2 = desired_intensity;

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

        % transmission setups are acoustic-only (no thermal properties in
        % tissue_setups.csv)
        if startsWith(setup_name, 'medium_trans_')
            parameters.modules.run_heating_sims = 0;
        end

        parameters.simulation.code_type = code_type;
        parameters.simulation.medium = 'phantom';
        parameters.grid.axisymmetric = 1;
        % 880 thermal steps, <= 3.8 GB peak RAM. CPU: 70-105+ min on dccn-c07x (node-load
        % dependent; 1:45 timed out 4/24 jobs on 2026-10-28). GPU: ~3 min.
        if contains(code_type, 'gpu')
            parameters.hpc.timelimit = '00:30:00';
        else
            parameters.hpc.timelimit = '03:00:00';
        end
        parameters.hpc.memorylimit = 8;

        % specify GPU setup if requested
        if contains(parameters.simulation.code_type, 'gpu')
            parameters.hpc.gpu = "gpu:1";
            parameters.hpc.partition = "gpu40g";
        end

        % Run the simulation
        parameters.subject_id = subject_id;
        parameters.platform = pn.platform;
        prestus_pipeline_start(parameters);

    end
    end
end
