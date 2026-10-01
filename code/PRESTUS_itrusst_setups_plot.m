% PRESTUS_ITRUSST_SETUPS_PLOT
%
% Parallel-coordinates plots of the brain and skull properties per setup:
% Setup 1 / Setup 2 (literature bounds; data/tissue_parameters.csv +
% tissue_setups.csv, via itrusst_tissue_properties.m) and the tool defaults
% PRESTUS / BabelBrain / k-Plan (data/configs/config_medium_<tool>.yaml).
% Attenuation is shown at the protocol frequency (500 kHz), i.e.
% alpha_coeff * f^alpha_power in dB/cm. k-Plan does not report perfusion
% (its config uses 0), so it is shown as missing.
%
% Per tissue, two figures: all axes min (bottom) to max (top), and a version
% with the axes flipped where Setup 1 uses the bottom bound, so that Setup 1
% and Setup 2 trace the upper and lower envelope.
%
% Output: figures/setups/itrusst_table_parallel_<tissue>[_flipped].{png,eps}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath); % PRESTUS yaml reader
pn.figures = fullfile(rootpath, 'figures', 'setups');
if ~isfolder(pn.figures); mkdir(pn.figures); end

%% setup definitions

freq_mhz = 0.5; % itrusst protocol 1

property_fields = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                   'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_labels = {'Sound Speed', 'Density', 'Att. (dB/cm)', 'Thermal conduct.', ...
                   'Specific Heat', 'Perf.', 'Abs. %'};

setup_labels = {'Setup 1', 'Setup 2', 'PRESTUS', 'BabelBrain', 'k-Plan'};
setup_colors = {[0.70 0.00 0.00], [0.00 0.60 0.50], [0.00 0.45 0.74], ...
                [0.49 0.18 0.56], [0.47 0.67 0.19]};

media = {itrusst_tissue_properties('setup1', pn.data), itrusst_tissue_properties('setup2', pn.data)};
tool_configs = {'PRESTUS', 'BabelBrain', 'kPlan'};
for i = 1:numel(tool_configs)
    cfg = yaml.loadFile(fullfile(pn.configs, ['config_medium_' tool_configs{i} '.yaml']));
    media{end+1} = cfg.medium_properties;
end

setups_table = readtable(fullfile(pn.data, 'tissue_setups.csv'), 'TextType', 'char');

%% plot per tissue

tissue_labels = {'brain', 'skull'};
for t = 1:numel(tissue_labels)
    tissue = tissue_labels{t};

    data = nan(numel(media), numel(property_fields));
    for s = 1:numel(media)
        m = media{s}.(tissue);
        for p = 1:numel(property_fields)
            data(s, p) = m.(property_fields{p});
        end
        data(s, strcmp(property_fields, 'alpha_coeff')) = m.alpha_coeff * freq_mhz^m.alpha_power;
    end
    data(strcmp(setup_labels, 'k-Plan'), strcmp(property_fields, 'perfusion')) = NaN; % not reported

    make_parallel_plot(data, property_labels, setup_labels, setup_colors, ...
        [], pn, ['itrusst_table_parallel_' tissue]);

    % flip the axes on which Setup 1 uses the bottom bound
    flip_cols = [];
    for p = 1:numel(property_fields)
        row = strcmp(setups_table.tissue, tissue) & strcmp(setups_table.property, property_fields{p});
        if strcmp(setups_table.setup1{row}, 'bottom')
            flip_cols(end+1) = p; %#ok<AGROW>
        end
    end
    make_parallel_plot(data, property_labels, setup_labels, setup_colors, ...
        flip_cols, pn, ['itrusst_table_parallel_' tissue '_flipped']);
end
