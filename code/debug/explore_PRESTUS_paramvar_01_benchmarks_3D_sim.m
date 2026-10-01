% PRESTUS_PARAMVAR_01_BENCHMARKS_3D_SIM
%
% Full 3-D (non-axisymmetric) benchmark simulations using the 3-D SimNIBS
% phantoms produced by explore_createPhantom3D.m (sub-101: water/tissue, sub-102: skull/brain).
%
% This is the 3-D counterpart of PRESTUS_PARAMVAR_01_BENCHMARKS_SIM, which
% uses axisymmetric 2-D grids.  Differences:
%   - Phantom folders:  m2m_sub-001_3D  /  m2m_sub-002_3D
%   - Config:           config_itrusst_protocol1_3D.yaml
%                         grid.axisymmetric = 0
%                         grid.default_dims = [140, 140, 280]
%                         transducer.trans_pos = [70, 70, 2]
%                         transducer.focus_pos = [70, 70, 130]
%   - Output folder:    data/tussim/itrusst_protocol1a_3D
%
% Run explore_createPhantom3D.m first if the phantom files do not yet exist.
%
% See also: CREATEPHANTOM3D, PRESTUS_PARAMVAR_01_BENCHMARKS_SIM

clear all; close all; clc;

user = 'julkos';

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, user);

%% Variant definitions (matching the 2-D study — baseline only for now)

variant_names     = {'itrusst_protocol1a_3D'};
variant_configs   = {'itrusst_protocol1_3D'};
variant_intensity = [60];

all_subjects = [101, 102];   % 3D phantoms: sub-101 = water/tissue, sub-102 = skull/brain

for subject_id = all_subjects

    for i_variant = 1:numel(variant_names)
        setup_name        = variant_configs{i_variant};
        variant_name      = variant_names{i_variant};
        desired_intensity = variant_intensity(i_variant);

        %% load parameters

        cd(fullfile(rootpath, 'data', 'configs'));
        parameters = load_parameters(['config_', setup_name, '.yaml']);

        if ismac
            parameters.startup.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
        else
            parameters.startup.simnibs_bin_path = pn.simnibs;
        end

        pn.path_sim = fullfile(rootpath, 'data', 'tussim', variant_name);
        if ~isfolder(pn.path_sim); mkdir(pn.path_sim); end

        parameters.hpc.ld_library_path         = "/opt/gcc/7.2.0/lib64";
        parameters.path.anat                   = pn.path_seg;
        parameters.path.seg                    = pn.path_seg;
        parameters.path.t1_pattern             = fullfile(sprintf('m2m_sub-%03d', subject_id), 'T1.nii.gz');
        parameters.path.t2_pattern             = fullfile(sprintf('m2m_sub-%03d', subject_id), 'T2_reg.nii.gz');
        parameters.path.sim                    = pn.path_sim;
        parameters.startup.paths_to_add        = {pn.kwave, pn.minimize};
        parameters.pct.enabled                 = 0;

        pn.outputs_folder = fullfile(pn.path_sim, sprintf('sub-%03d', subject_id));
        if ~isfolder(pn.outputs_folder); mkdir(pn.outputs_folder); end

        %% run phantom simulation — mid-range tissue properties (baseline only)

        parameters.simulation.medium      = 'phantom';
        parameters.simulation.interactive = 0;
        parameters.io.overwrite_files     = 'never';
        parameters.io.overwrite_simnibs   = 0;
        parameters.io.save_heatingvideo   = 0;
        parameters.io.save_matrices       = 0;
        parameters.io.output_affix        = '';

        parameters.simulation.code_type   = 'cpp_gpu';
        parameters.hpc.gpu                = "gpu:1";
        parameters.hpc.partition          = "gpu40g";
        parameters.simulation.precision   = 'single';
        parameters.subject_id             = subject_id;

        prestus_pipeline_start(parameters);
    end
end
