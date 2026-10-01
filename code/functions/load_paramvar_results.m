function [results, xlabels] = load_paramvar_results(data_dir, subject_id, ...
    tissue_labels, tissue_suffixes, property_labels, property_suffixes, sweep_values, property_symbols)
% Load parameter variation simulation results from CSV files.
%
% Inputs:
%   data_dir          - folder containing sub-XXX_phantom*.csv files
%   subject_id        - numeric subject ID (e.g. 1, 2, 303)
%   tissue_labels     - cell array of tissue names, e.g. {'brain','skull'}
%   tissue_suffixes   - cell array of file suffixes,  e.g. {'_b','_s'}
%   property_labels   - cell array of property names
%   property_suffixes - cell array of file suffixes for each property
%   sweep_values      - numeric vector of sweep level indices (e.g. 1:2)
%   property_symbols  - (optional) cell array of TeX symbols for each property;
%                       if provided, labels are formatted as 'sym_{tissue}'
%
% Outputs:
%   results  - struct with fields: sweep, Ipa_target, maxT (rise from 37 C),
%              Isppa_brain, riseT_brain, Isppa_skull, riseT_skull,
%              real_focal_distance_mm, focal_vol_3dB_mm3 (ITRUSST −3 dB /
%              half-intensity focal volume, i.e. the volumetric FWHM)
%              — each [n_sweep x n_conds]
%   xlabels  - cell array of condition labels (length n_conds)

if nargin < 8
    property_symbols = [];
end

n_conds = length(tissue_labels) * length(property_labels);
n_sweep = length(sweep_values);
sub_str = sprintf('sub-%03d', subject_id);

metric_cols = {'Ipa_target', 'maxT', 'Isppa_brain', 'riseT_brain', ...
               'Isppa_skull', 'riseT_skull', 'real_focal_distance_mm', ...
               'focal_vol_3dB_mm3'};

results  = struct();
xlabels  = cell(1, n_conds);
for m = 1:length(metric_cols)
    results.(metric_cols{m}) = nan(n_sweep, n_conds);
end
results.sweep = nan(n_sweep, n_conds);

i_cond = 0;
for t = 1:length(tissue_labels)
    tissue_suffix = tissue_suffixes{t};

    for p = 1:length(property_labels)
        i_cond = i_cond + 1;
        prop_suffix = property_suffixes{p};
        tissue_sup = tissue_suffix(2:end);
        if ~isempty(property_symbols)
            xlabels{i_cond} = [property_symbols{p} '_{' tissue_sup '}'];
        else
            xlabels{i_cond} = [tissue_sup '_' property_labels{p}];
        end
        results.sweep(:, i_cond) = sweep_values(:);

        for v = 1:length(sweep_values)
            affix    = [tissue_suffix prop_suffix num2str(sweep_values(v))];
            filename = fullfile(data_dir, [sub_str '_phantom' affix '.csv']);

            if ~isfile(filename)
                warning('File not found: %s', filename);
                continue;
            end

            T    = readtable(filename);
            vars = T.Properties.VariableNames;

            for m = 1:length(metric_cols)
                col = metric_cols{m};
                if ismember(col, vars)
                    results.(col)(v, i_cond) = T.(col)(1);
                else
                    warning('File %s missing %s.', filename, col);
                end
            end
        end
    end
end

results.maxT = results.maxT - 37;
