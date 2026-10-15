function p = itrusst_input_path(default_path, env_var, what)
% Path of an input that is not tracked in the repository (collaborator
% results, external literature tables). The environment variable env_var
% overrides the default location; an error names both if neither exists.
%
% Inputs:
%   default_path - expected location (file or folder) in the project tree
%   env_var      - environment variable with a user-provided location,
%                  e.g. setenv('ITRUSST_BABELBRAIN_DIR', '/path/to/csvs')
%   what         - short description used in the error message
%
% Output:
%   p            - existing file or folder

p = getenv(env_var);
source = sprintf('environment variable %s', env_var);
if isempty(p)
    p = default_path;
    source = 'default location';
end

if ~isfile(p) && ~isfolder(p)
    error('itrusst:missingInput', ['%s not found at %s (%s).\n' ...
        'This input is not tracked in the repository. Provide its location with\n' ...
        '  setenv(''%s'', ''<path>'')\n' ...
        'before running the script.'], what, p, source, env_var);
end
end
