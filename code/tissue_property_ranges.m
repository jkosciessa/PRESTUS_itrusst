%% Tissue Property Ranges Plot with plot_metric_with_patches

% 1. Load & process data (from your working results table)
pn.data = '/Volumes/2425122.01/PRESTUS_itrusst/data';
params_table = readtable(fullfile(pn.data, 'tissue_parameters.csv'));

tissue_labels = {'brain', 'skull'};
property_labels = {'sound_speed', 'density', 'alpha_0_true', 'thermal_conductivity', ...
                   'specific_heat_capacity', 'perfusion', 'absorption_fraction'};

property_names = {'sound_speed', 'density', 'attenuation', 'thermal_conductivity', ...
                   'heat_capacity', 'perfusion', 'absorption'};

% Build results table (your fixed version)
results = table();
for t = 1:length(tissue_labels)
    tissue = tissue_labels{t};
    tissue_mask = strcmp(params_table.tissue, tissue);
    for p = 1:length(property_labels)
        prop = property_labels{p};
        prop_name = property_names{p};
        mid_val = params_table.(prop)(tissue_mask & strcmp(params_table.level, 'mid'));
        top_val = params_table.(prop)(tissue_mask & strcmp(params_table.level, 'top'));
        bot_val = params_table.(prop)(tissue_mask & strcmp(params_table.level, 'bottom'));
        
        if length(mid_val) == 1 && length(top_val) == 1 && length(bot_val) == 1
            pct_top = 100 * (top_val - mid_val) / mid_val;
            pct_bot = 100 * (bot_val - mid_val) / mid_val;
            
            new_row = {tissue, prop_name, mid_val, top_val, bot_val, pct_top, pct_bot};
            temp_table = cell2table(new_row, 'VariableNames', ...
                {'Tissue','Property','MID','TOP','BOTTOM','PctTopMid','PctBotMid'});
            results = [results; temp_table];
        end
    end
end

% 2. Prepare for plot_metric_with_patches (14 properties x 2 levels)
data_matrix = [results.PctTopMid, results.PctBotMid]';  % [14 x 2]

% X-labels: Tissue_Property format
xlabels = strcat(results.Tissue, '_', regexprep(results.Property, '_', '-'));

% 3. Create output dir
pn.figures = fullfile(pn.data, 'figures');
if ~exist(pn.figures, 'dir'), mkdir(pn.figures); end

% 4. Plot ALL 14 on single axis
plot_metric_with_patches(double(data_matrix), xlabels, '% Change vs MID', ...
    'all_tissue_ranges.png', [], [1; 2], ...
    {[0 0 0], [1 1 1]}, pn);