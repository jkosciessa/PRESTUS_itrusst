% PRESTUS_PARAMVAR_01_BENCHMARKS_025MM_SIM
%
% 2-D axisymmetric benchmark simulations at 0.5 mm resolution (6 PPW at
% 500 kHz) using phantoms produced by explore_createPhantom025mm.m.
%
% Phantom folders:  m2m_sub-201  /  m2m_sub-202
% Config:           config_itrusst_protocol1_025mm.yaml
%                     grid.resolution_mm = 0.5
%                     grid.default_dims  = [140, 240]
%                     grid.pml_size      = 40
%                     transducer.trans_pos = [70, 1]
% Output folder:    data/tussim/itrusst_protocol1a_025mm
%
% Run explore_createPhantom025mm.m first if the phantom files do not yet exist.
%
% See also: CREATEPHANTOM025MM, PRESTUS_PARAMVAR_01_BENCHMARKS_SIM,
%           PRESTUS_PARAMVAR_01_BENCHMARKS_3D_SIM

clear all; close all; clc;

user = 'julkos';

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, user);

%% Variant definitions

variant_names     = {'itrusst_protocol1a_05mm'};
variant_configs   = {'itrusst_protocol1_05mm'};
variant_intensity = [60];

all_subjects = 002;%[001, 002];   % 0.5 mm phantoms: sub-201 = water/tissue, sub-202 = skull/brain

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

        %% run phantom simulation

        parameters.simulation.medium      = 'phantom';
        parameters.simulation.interactive = 0;
        parameters.io.overwrite_files     = 'always';
        parameters.io.overwrite_simnibs   = 0;
        parameters.io.save_heatingvideo   = 0;
        parameters.io.save_matrices       = 0;
        parameters.io.output_affix        = '';

        parameters.simulation.code_type   = 'cpp_gpu';
        parameters.simulation.precision   = 'single';
        parameters.subject_id             = subject_id;
        parameters.platform = 'matlab';

        % parameters.grid.source_cfl = 0.01;
        parameters.grid.pml_size = 'auto';

        prestus_pipeline_start(parameters);
    end
end
