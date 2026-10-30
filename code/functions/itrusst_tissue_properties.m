function medium = itrusst_tissue_properties(setup_name, data_dir)
% Tissue properties of a literature setup, read from the project CSVs.
%
%   data/tissue_parameters.csv - bottom/mid/top value per tissue x property
%                                (alpha_coeff in dB/(MHz^y cm), y = alpha_power;
%                                brain absorption as absorption_fraction, skull
%                                and water absorption as absorption_coeff;
%                                water has a single 'mid' row)
%   data/tissue_setups.csv     - which level (bottom/mid/top) each setup uses per
%                                tissue x property; one column per setup
%                                (heat_max, heat_min, trans_max, trans_min);
%                                an empty entry = not specified (trans_*:
%                                acoustic-only setups, thermal properties
%                                not simulated) and is filled with 'mid' so
%                                that the medium is complete
%
% Inputs:
%   setup_name - a setup column of tissue_setups.csv (level per property),
%                'reference' (= mid), or a level name ('bottom', 'mid', 'top')
%   data_dir   - folder containing the CSVs (pn.data)
%
% Output:
%   medium     - struct medium.(water|brain|skull).(property), as used in
%                parameters.medium_properties; water is always the 'mid' row;
%                absorption resolved to absorption_coeff (itrusst_resolve_absorption)

params_table = readtable(fullfile(data_dir, 'tissue_parameters.csv'), 'TextType', 'char');
setups_table = readtable(fullfile(data_dir, 'tissue_setups.csv'), 'TextType', 'char');

property_labels = setdiff(params_table.Properties.VariableNames, {'tissue', 'level'}, 'stable');
tissue_labels   = unique(params_table.tissue, 'stable');
setup_columns   = setdiff(setups_table.Properties.VariableNames, {'tissue', 'property'}, 'stable');

if strcmp(setup_name, 'reference')
    setup_name = 'mid';
end
if ~ismember(setup_name, [setup_columns, {'bottom', 'mid', 'top'}])
    error('Unknown setup ''%s''; tissue_setups.csv defines %s.', setup_name, strjoin(setup_columns, ', '));
end

medium = struct();
for t = 1:numel(tissue_labels)
    tissue = tissue_labels{t};
    for p = 1:numel(property_labels)
        property = property_labels{p};
        if strcmp(tissue, 'water')
            level = 'mid';
        elseif ismember(setup_name, setup_columns)
            row = strcmp(setups_table.tissue, tissue) & strcmp(setups_table.property, property);
            if any(row)
                level = setups_table.(setup_name){row};
                if isempty(level); level = 'mid'; end % not specified
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

% brain: absorption fraction x the attenuation of this setup;
% skull/water: absorption coefficient, capped at the attenuation
medium = itrusst_resolve_absorption(medium);
