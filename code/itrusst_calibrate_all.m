clear all; close all; clc;

user = 'julkos'

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, user);

%% define variants
% Calibration is transducer-specific and does not depend on subject anatomy.
variant_configs   = {'itrusst_protocol'};
variant_intensity = [60];

% Use a fixed subject_id for the reference free-water simulation
subject_id = 1;

for i_variant = 1:length(variant_configs)
    setup_name        = variant_configs{i_variant};
    desired_intensity = variant_intensity(i_variant);

    cd(fullfile(rootpath, 'data', 'configs'))
    parameters = load_parameters(['config_', setup_name, '.yaml']);

    if ismac
        parameters.startup.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
    else
        parameters.startup.simnibs_bin_path = pn.simnibs;
    end
    pn.path_sim = fullfile(rootpath, 'data', 'tussim', setup_name);
        if ~exist(pn.path_sim); mkdir(pn.path_sim); end
    parameters.hpc.ld_library_path = "/opt/gcc/7.2.0/lib64";
    parameters.path.anat = pn.path_seg;
    parameters.path.seg  = pn.path_seg;
    parameters.path.t1_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
    parameters.path.t2_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
    parameters.path.sim = pn.path_sim;
    parameters.startup.paths_to_add = {pn.kwave, pn.minimize};
    pn.outputs_folder = fullfile(pn.path_sim, sprintf('sub-%03d', subject_id));
        if ~exist(pn.outputs_folder); mkdir(pn.outputs_folder); end
    parameters.pct.enabled = 0;
    parameters.platform = 'slurm';

    itrusst_calibrate(pn, parameters, subject_id, desired_intensity, setup_name);

end
