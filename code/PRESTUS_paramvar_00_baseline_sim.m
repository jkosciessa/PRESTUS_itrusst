clear all; close all; clc;

user = 'julkos';

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
rootpath = fileparts(pathstr);

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, user);

%% settings

all_subjects      = [001, 002, 352];   % as PRESTUS_paramvar_01_benchmarks_sim
setup_name        = 'itrusst_protocol';
variant_name      = 'itrusst_protocol1a';

for subject_id = all_subjects

    % skip subjects whose baseline already has thermal output
    baseline_csv = fullfile(rootpath, 'data', 'tussim', variant_name, sprintf('sub-%03d', subject_id), ...
        sprintf('sub-%03d_phantom_baseline.csv', subject_id));
    if isfile(baseline_csv) && ismember('maxT', readtable(baseline_csv).Properties.VariableNames)
        fprintf('sub-%03d: baseline complete, skipped.\n', subject_id);
        continue;
    end

    %% load parameter file and adjust paths

    cd(fullfile(rootpath, 'data', 'configs'));
    parameters = load_parameters(['config_', setup_name, '.yaml']);


    pn.path_sim = fullfile(rootpath, 'data', 'tussim', variant_name);
    if ~exist(pn.path_sim, 'dir'); mkdir(pn.path_sim); end

    parameters.hpc.ld_library_path = "/opt/gcc/7.2.0/lib64";
    parameters.path.anat            = pn.path_seg;
    parameters.path.seg             = pn.path_seg;
    parameters.path.t1_pattern      = fullfile(sprintf('m2m_sub-%03d', subject_id), 'T1.nii.gz');
    parameters.path.t2_pattern      = fullfile(sprintf('m2m_sub-%03d', subject_id), 'T2_reg.nii.gz');
    parameters.path.sim             = pn.path_sim;
    parameters.startup.paths_to_add = {pn.kwave, pn.minimize};

    pn.outputs_folder = fullfile(pn.path_sim, sprintf('sub-%03d', subject_id));
    if ~exist(pn.outputs_folder, 'dir'); mkdir(pn.outputs_folder); end

    %% apply mid-level (default) tissue properties
    % tissue_parameters.csv level='mid' — the baseline of the paramvar sweep;
    % the config holds the same values, set here from the CSV so they stay in sync.

    tissue_labels   = {'brain', 'skull'};
    property_labels = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                       'specific_heat_capacity', 'perfusion', 'absorption_coeff', ...
                       'absorption_fraction'}; % brain absorption as fraction, skull as coefficient

    params_table = readtable(fullfile(pn.data, 'tissue_parameters.csv'));
    for t = 1:length(tissue_labels)
        for p = 1:length(property_labels)
            parameters.medium_properties.(tissue_labels{t}).(property_labels{p}) = ...
                params_table.(property_labels{p})(strcmp(params_table.tissue, tissue_labels{t}) & ...
                                                   strcmp(params_table.level, 'mid'));
        end
    end
    parameters.medium_properties = itrusst_resolve_absorption(parameters.medium_properties);

    %% simulation settings

    parameters.simulation.medium      = 'phantom';
    parameters.simulation.interactive = 0;
    parameters.simulation.code_type   = 'cpp_cpu'; % as the paramvar_01 sweep it normalises
    parameters.simulation.precision   = 'single';

    parameters.io.overwrite_files    = 'never';
    parameters.io.overwrite_simnibs  = 0;
    parameters.io.save_heatingvideo  = 0;
    parameters.io.save_matrices      = 0;

    parameters.grid.resolution_mm  = 0.5;
    parameters.grid.default_dims   = [70, 120] / parameters.grid.resolution_mm;
    parameters.grid.pml_size       = [20]  / parameters.grid.resolution_mm; % 20 mm PML, matches primary pipeline configs
    parameters.grid.axisymmetric   = 1;

    parameters.pct.enabled = 0;

    parameters.platform        = pn.platform; % scheduler if available, else MATLAB (setup_itrusst_paths)
    parameters.hpc.timelimit   = '03:00:00'; % as PRESTUS_paramvar_01_benchmarks_sim
    parameters.hpc.memorylimit = 8;          % peak RAM <= 3.8 GB

    % all-mid baseline (counterpart to the bottom/top sweep of paramvar_01)
    parameters.io.output_affix = '_baseline';

    parameters.subject_id = subject_id;

    %% run

    prestus_pipeline_start(parameters);

end
