function [results, ref, xlabels] = load_babelbrain_benchmark2(data_dir, ppw, ...
    tissue_labels, property_labels, property_symbols)
% BabelBrain results of the one-at-a-time parameter variation on benchmark 2
% (skull + brain phantom, = PRESTUS sub-002), provided by Sam Pichardo as
% data/pichardo/Benchmark2-<Mon-dd-yyyy>-<ppw>PPW.csv (newest date is used).
%
% Same layout as load_paramvar_results: conditions are tissue x property
% (tissue outer loop), rows are the LOW / HIGH level (= sweep level 1 / 2,
% bottom / top of data/tissue_parameters.csv).
%
% Inputs:
%   data_dir          - folder with the BabelBrain CSVs (data/pichardo)
%   ppw               - points per wavelength of the run (6 or 9)
%   tissue_labels     - e.g. {'brain','skull'}
%   property_labels   - PRESTUS property names: sound_speed, density,
%                       alpha_coeff, thermal_conductivity,
%                       specific_heat_capacity, perfusion, absorption_fraction
%   property_symbols  - (optional) TeX symbols, labels become 'sym_{b|s}'
%
% Outputs:
%   results  - struct, each field [2 x n_conds]:
%                Isppa_brain   Max. Intensity Brain (W/cm^2)
%                MI_brain      MI
%                riseT_brain   Max. Temperature Brain (rise, °C)
%                riseT_skull   Max. Temperature Skull (rise, °C)
%                riseT_target  Max. Temperature Target (rise, °C)
%                CEM43_brain   Thermal Dose Brain
%                CEM43_skull   Thermal Dose Skull
%   ref      - same fields as scalars, from the 'Reference' row (all 'mid')
%   xlabels  - condition labels (length n_conds)

if nargin < 5
    property_symbols = [];
end

csv_file = find_csv(data_dir, ppw);

% BabelBrain column per output field
cols = {'Isppa_brain',  'Max. Intensity Brain'
        'MI_brain',     'MI'
        'riseT_brain',  'Max. Temperature Brain'
        'riseT_skull',  'Max. Temperature Skull'
        'riseT_target', 'Max. Temperature Target'
        'CEM43_brain',  'Thermal Dose Brain'
        'CEM43_skull',  'Thermal Dose Skull'};

% BabelBrain name per PRESTUS property
bb_prop = struct('sound_speed', 'SoS', 'density', 'Density', 'alpha_coeff', 'Attenuation', ...
    'thermal_conductivity', 'Thermal Conductivity', 'specific_heat_capacity', 'Specific Heat', ...
    'perfusion', 'Perfusion', 'absorption_fraction', 'Absorption');
levels = {'LOW', 'HIGH'};

T = readtable(csv_file, 'VariableNamingRule', 'preserve', 'TextType', 'char', 'Delimiter', ',');
type = strtrim(T.Type);

n_conds = numel(tissue_labels) * numel(property_labels);
results = struct();
for c = 1:size(cols, 1)
    results.(cols{c, 1}) = nan(numel(levels), n_conds);
end
xlabels = cell(1, n_conds);

i_cond = 0;
for t = 1:numel(tissue_labels)
    tissue = tissue_labels{t};
    for p = 1:numel(property_labels)
        i_cond = i_cond + 1;
        if ~isempty(property_symbols)
            xlabels{i_cond} = [property_symbols{p} '_{' tissue(1) '}'];
        else
            xlabels{i_cond} = [tissue(1) '_' property_labels{p}];
        end
        for v = 1:numel(levels)
            row_name = sprintf('%s%s %s %s', upper(tissue(1)), tissue(2:end), ...
                bb_prop.(property_labels{p}), levels{v});
            row = find(strcmp(type, row_name), 1);
            if isempty(row)
                warning('%s: no row ''%s''.', csv_file, row_name);
                continue;
            end
            for c = 1:size(cols, 1)
                results.(cols{c, 1})(v, i_cond) = T.(cols{c, 2})(row);
            end
        end
    end
end

ref = struct();
row = find(strcmp(type, 'Reference'), 1);
for c = 1:size(cols, 1)
    ref.(cols{c, 1}) = NaN;
    if ~isempty(row); ref.(cols{c, 1}) = T.(cols{c, 2})(row); end
end
if isempty(row); warning('%s: no ''Reference'' row.', csv_file); end
end

function csv_file = find_csv(data_dir, ppw)
% newest Benchmark2-<Mon-dd-yyyy>-<ppw>PPW.csv
files = dir(fullfile(data_dir, sprintf('Benchmark2-*-%dPPW.csv', ppw)));
if isempty(files)
    error('itrusst:missingInput', ['No BabelBrain benchmark 2 results for %d PPW ' ...
        '(Benchmark2-<Mon-dd-yyyy>-%dPPW.csv) in %s.\nThese results are not tracked ' ...
        'in the repository; setenv(''ITRUSST_BABELBRAIN_DIR'', ''<folder>'') points the ' ...
        'scripts to another copy.'], ppw, ppw, data_dir);
end
dates = NaT(numel(files), 1);
for f = 1:numel(files)
    tok = regexp(files(f).name, '^Benchmark2-(\w{3}-\d{1,2}-\d{4})-', 'tokens', 'once');
    if ~isempty(tok)
        dates(f) = datetime(tok{1}, 'InputFormat', 'MMM-dd-yyyy', 'Locale', 'en_US');
    end
end
[~, i_new] = max(dates);
if isempty(i_new) || isnat(dates(i_new)); i_new = numel(files); end
csv_file = fullfile(data_dir, files(i_new).name);
end
