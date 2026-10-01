function pn = setup_itrusst_paths(rootpath, user)
% Initialize PRESTUS toolbox paths and add them to the MATLAB path.
%
% Inputs:
%   rootpath - root directory of the PRESTUS_itrusst project
%   user     - (optional) username string; sets pn.simnibs for known users

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

if nargin >= 2 && strcmp(user, 'julkos')
    pn.simnibs = fullfile('/home', 'neuromod', 'julkos', '.conda', 'envs', 'simnibs_env', 'bin');
end
