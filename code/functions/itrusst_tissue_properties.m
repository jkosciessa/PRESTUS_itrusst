function medium = itrusst_tissue_properties(setup_name, data_dir)
% Tissue properties of a literature setup, read from the project CSVs.
%
%   data/tissue_parameters.csv - bottom/mid/top value per tissue x property
%                                (alpha_coeff in dB/(MHz^y cm), y = alpha_power;
%                                water has a single 'mid' row)
%   data/tissue_setups.csv     - which bound (bottom/top) setup1/setup2 use per
%                                tissue x property
%
% Inputs:
%   setup_name - 'setup1' or 'setup2' (bound per property from tissue_setups.csv),
%                'reference' (= mid), or a level name ('bottom', 'mid', 'top')
%   data_dir   - folder containing the CSVs (pn.data)
%
% Output:
%   medium     - struct medium.(water|brain|skull).(property), as used in
%                parameters.medium_properties; water is always the 'mid' row

params_table = readtable(fullfile(data_dir, 'tissue_parameters.csv'), 'TextType', 'char');
setups_table = readtable(fullfile(data_dir, 'tissue_setups.csv'), 'TextType', 'char');

property_labels = setdiff(params_table.Properties.VariableNames, {'tissue', 'level'}, 'stable');
tissue_labels   = unique(params_table.tissue, 'stable');

if strcmp(setup_name, 'reference')
    setup_name = 'mid';
end

medium = struct();
for t = 1:numel(tissue_labels)
    tissue = tissue_labels{t};
    for p = 1:numel(property_labels)
        property = property_labels{p};
        if strcmp(tissue, 'water')
            level = 'mid';
        elseif ismember(setup_name, {'setup1', 'setup2'})
            row = strcmp(setups_table.tissue, tissue) & strcmp(setups_table.property, property);
            if any(row)
                level = setups_table.(setup_name){row};
            else
                level = 'mid'; % e.g. alpha_power: not varied between setups
            end
        else
            level = setup_name;
        end
        value = params_table.(property)(strcmp(params_table.tissue, tissue) & strcmp(params_table.level, level));
        if numel(value) ~= 1
            error('No unique %s %s value at level ''%s'' in tissue_parameters.csv.', tissue, property, level);
        end
        medium.(tissue).(property) = value;
    end
end
