function pn = setup_itrusst_paths(rootpath, user)
% Initialize PRESTUS toolbox paths and add them to the MATLAB path.
%
% Inputs:
%   rootpath - root directory of the PRESTUS_itrusst project
%   user     - (optional) username string; sets pn.simnibs for known users
%
% Also sets pn.platform (job submission) and pn.simnibs ('' unless user is known).

% On macOS the volume may be accessed via /project/ alias; normalize to /Volumes/
if ismac
    rootpath = strrep(rootpath, '/project/', '/Volumes/');
end

pn.prestus         = fullfile(rootpath, 'tools', 'PRESTUS'); % submodule, pulled from ../v05/tools/PRESTUS
pn.prestus_fun     = fullfile(pn.prestus, 'functions');
pn.prestus_tools   = fullfile(pn.prestus, 'external');
pn.kwave           = fullfile(pn.prestus_tools, 'k-wave', 'k-Wave');
pn.minimize        = fullfile(pn.prestus_tools, 'FEX-minimize');
pn.nifti           = fullfile(rootpath, '..', 'v05', 'tools', 'nifti_toolbox');
pn.configs         = fullfile(rootpath, 'data', 'configs');
pn.data_path       = fullfile(rootpath, 'data', 'bids');
pn.path_seg        = fullfile(rootpath, 'data', 'simnibs');
pn.data            = fullfile(rootpath, 'data');
pn.code            = fullfile(rootpath, 'code');

addpath(genpath(pn.prestus));
addpath(genpath(pn.prestus_fun));
addpath(genpath(pn.prestus_tools));
addpath(pn.kwave);
addpath(pn.minimize);
addpath(pn.nifti);
addpath(fullfile(pn.code, 'functions'));

% the configs set thermal.absorption_model: coefficient; an older PRESTUS ignores
% that field and deposits 100% of attenuation as heat
if exist('thermal_absorption_model', 'file') ~= 2
    error('setup_itrusst_paths:noAbsorptionCoeff', ...
        ['PRESTUS in %s does not support thermal.absorption_model (functions/thermal/' ...
         'thermal_absorption_model.m missing); update the submodule from ../v05/tools/PRESTUS.'], pn.prestus);
end

% SimNIBS binaries: only needed to segment new heads, not for the phantoms here
pn.simnibs = '';
if nargin >= 2 && strcmp(user, 'julkos')
    if ismac
        pn.simnibs = '/Users/julian.kosciessa/SimNIBS/bin/';
    else
        pn.simnibs = fullfile('/home', 'neuromod', 'julkos', '.conda', 'envs', 'simnibs_env', 'bin');
    end
end

% job submission: the scheduler on the PATH (slurm/qsub), otherwise in MATLAB
% (PRESTUS platform 'auto' errors without a scheduler);
% setenv('ITRUSST_PLATFORM', 'matlab') runs locally on a cluster node
pn.platform = getenv('ITRUSST_PLATFORM');
if isempty(pn.platform)
    try
        pn.platform = hpc_detect_system();
    catch
        pn.platform = 'matlab';
    end
end
fprintf('Job submission platform: %s\n', pn.platform);
