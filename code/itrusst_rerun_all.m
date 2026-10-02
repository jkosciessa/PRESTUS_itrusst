function itrusst_rerun_all(stage)
% Rerun the main analysis chain after a change of tissue_parameters.csv.
%
%   itrusst_rerun_all('prepare')    % local, minutes: medium configs, thickness
%                                   % phantoms, tissue-range and setup figures
%   itrusst_rerun_all('calibrate')  % cluster, blocks until done: free-water
%                                   % calibration, copied into the protocol configs
%   itrusst_rerun_all('simulate')   % cluster: submits all simulations and returns
%   itrusst_rerun_all('plot')       % after all simulation jobs have finished
%
% Run the stages in this order. Before the first stage, move the previous
% outputs aside with itrusst_backup_outputs(rootpath, '<yymmdd>'), as several
% simulations use overwrite_files = 'never' and would skip existing outputs.
%
% Each step is one of the scripts in code/, run unchanged (they clear the
% workspace, so each runs in its own function workspace). A failing step is
% reported and the remaining steps of the stage still run.

if nargin < 1
    error('Specify a stage: ''prepare'', ''calibrate'', ''simulate'' or ''plot''.');
end

user = 'julkos';

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
rootpath = fileparts(pathstr);

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, user);

switch stage
    case 'prepare'
        itrusst_write_medium_configs(pn);
        steps = {'PRESTUS_thickness_01_createPhantoms', ...
                 'tissue_property_ranges', ...
                 'PRESTUS_itrusst_setups_plot'};

    case 'calibrate'
        steps = {'itrusst_calibrate_all'};

    case 'simulate'
        cal_dir = fullfile(pn.data, 'calibration', 'PRESTUS_virtual_parameters');
        if isempty(dir(fullfile(cal_dir, 'itrusst_protocol-F*mm-I60wpercm2.yaml')))
            error('No calibration in %s: run itrusst_rerun_all(''calibrate'') first.', cal_dir);
        end
        steps = {'PRESTUS_paramvar_01_benchmarks_sim', ...
                 'PRESTUS_pipeline_01_benchmark1_sim', ...
                 'PRESTUS_thickness_02_sim', ...
                 'PRESTUS_thickness_04_density_sim', ...
                 'PRESTUS_paramvar_00_baseline_sim'}; % runs locally on the GPU (platform 'matlab')

    case 'plot'
        steps = {'PRESTUS_paramvar_02_benchmark1_plot', ...
                 'PRESTUS_paramvar_02_benchmark2_plot', ...
                 'PRESTUS_paramvar_02_benchmark_allvariants_plot', ...
                 'PRESTUS_paramvar_02_thermal_timeseries_plot', ...
                 'PRESTUS_paramvar_04_axial_profiles_2D', ...
                 'PRESTUS_itrusst_overview_profiles', ...
                 'PRESTUS_pipeline_02_benchmark1_plot', ...
                 'PRESTUS_thickness_03_plot', ...
                 'PRESTUS_thickness_05_density_plot'};

    otherwise
        error('Unknown stage ''%s''.', stage);
end

failed = {};
for i = 1:numel(steps)
    fprintf('\n=== [%s %d/%d] %s ===\n', stage, i, numel(steps), steps{i});
    try
        run_step(fullfile(rootpath, 'code', [steps{i} '.m']));
    catch err
        fprintf(2, 'FAILED %s: %s\n', steps{i}, err.message);
        failed{end+1} = steps{i}; %#ok<AGROW>
    end
    close all;
end

% the calibrated drive goes into the configs the simulations load
if strcmp(stage, 'calibrate') && isempty(failed)
    itrusst_apply_calibration(pn, 'itrusst_protocol', 60);
end

cd(rootpath);
if isempty(failed)
    fprintf('\nStage ''%s'' done (%d steps).\n', stage, numel(steps));
else
    fprintf(2, '\nStage ''%s'': %d of %d steps failed: %s\n', stage, numel(failed), numel(steps), strjoin(failed, ', '));
end
if strcmp(stage, 'simulate')
    fprintf('Run itrusst_rerun_all(''plot'') once all submitted jobs have finished.\n');
end

end

function run_step(script_file)
% own workspace: the scripts start with 'clear all'
run(script_file);
end
