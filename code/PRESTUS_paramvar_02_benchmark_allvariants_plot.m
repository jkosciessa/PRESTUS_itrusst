% Plot parameter variation results for all three protocol variants and all subjects.
% Variant data paths follow the naming from PRESTUS_paramvar_01_benchmarks_sim.m:
%   itrusst_protocol1a  (500 kHz, 60 W/cm2)
%   itrusst_protocol1b  (500 kHz, 30 W/cm2)
%   itrusst_protocol1c  (250 kHz, 60 W/cm2)
% Figures are saved in figures/<variant_name>/.

% get root path (script must be run from file location)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.code = fullfile(rootpath, 'code'); addpath(pn.code)

%% Configuration

variant_names    = {'itrusst_protocol1a', 'itrusst_protocol1b', 'itrusst_protocol1c'};
all_subjects     = [001, 002, 303];

% Tissue definitions (same for all subjects; sub-001 has brain only)
tissue_labels_sub001 = {'brain'};
tissue_suffixes_sub001 = {'_b'};

tissue_labels_sub002303 = {'brain', 'skull'};
tissue_suffixes_sub002303 = {'_b', '_s'};

% Property definitions (consistent with simulation script)
property_labels   = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                     'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};

sweep_values  = 1:2;
scatter_levels = [1, 2];
scatter_colors = {[1 1 1], [0 0 0]};

%% Loop over variants

for i_variant = 1:length(variant_names)
    variant = variant_names{i_variant};

    pn.figures = fullfile(rootpath, 'figures', variant);
    if ~exist(pn.figures, 'dir'); mkdir(pn.figures); endo

    %% Loop over subjects

    for subject_id = all_subjects
        sub_str = sprintf('sub-%03d', subject_id);
        pn.data = fullfile(rootpath, 'data', 'tussim', variant, sub_str);

        if subject_id == 1
            tissue_labels   = tissue_labels_sub001;
            tissue_suffixes = tissue_suffixes_sub001;
        else
            tissue_labels   = tissue_labels_sub002303;
            tissue_suffixes = tissue_suffixes_sub002303;
        end

        n_conds = length(tissue_labels) * length(property_labels);

        % Preallocate
        results       = struct();
        xlabels       = {};
        Isppa_brain   = nan(length(sweep_values), n_conds);
        riseT_brain   = nan(length(sweep_values), n_conds);
        real_focal_mm = nan(length(sweep_values), n_conds);

        has_skull = length(tissue_labels) > 1;
        if has_skull
            Isppa_skull = nan(length(sweep_values), n_conds);
            riseT_skull = nan(length(sweep_values), n_conds);
        end

        i_cond = 0;
        for t = 1:length(tissue_labels)
            tissue        = tissue_labels{t};
            tissue_suffix = tissue_suffixes{t};

            for p = 1:length(property_labels)
                i_cond      = i_cond + 1;
                property    = property_labels{p};
                prop_suffix = property_suffixes{p};

                isppa_vec       = nan(length(sweep_values), 1);
                maxT_vec        = nan(length(sweep_values), 1);
                maxIsppa_b_vec  = nan(length(sweep_values), 1);
                riseT_b_vec     = nan(length(sweep_values), 1);
                real_focal_vec  = nan(length(sweep_values), 1);
                if has_skull
                    maxIsppa_s_vec = nan(length(sweep_values), 1);
                    riseT_s_vec    = nan(length(sweep_values), 1);
                end

                for v = 1:length(sweep_values)
                    affix    = [tissue_suffix prop_suffix num2str(sweep_values(v))];
                    xlabels{i_cond} = [tissue_suffix(2:end) '_' property];
                    filename = fullfile(pn.data, [sub_str '_phantom_output_table' affix '.csv']);

                    if isfile(filename)
                        T = readtable(filename);

                        if ismember('Ipa_target', T.Properties.VariableNames)
                            isppa_vec(v) = T.Ipa_target(1);
                        else
                            warning('File %s missing Ipa_target.', filename);
                        end

                        if ismember('Isppa_brain', T.Properties.VariableNames)
                            maxIsppa_b_vec(v) = T.Isppa_brain(1);
                        else
                            warning('File %s missing Isppa_brain.', filename);
                        end

                        if ismember('maxT', T.Properties.VariableNames)
                            maxT_vec(v) = T.maxT(1);
                        else
                            warning('File %s missing maxT.', filename);
                        end

                        if ismember('riseT_brain', T.Properties.VariableNames)
                            riseT_b_vec(v) = T.riseT_brain(1);
                        else
                            warning('File %s missing riseT_brain.', filename);
                        end

                        if ismember('real_focal_distance_mm', T.Properties.VariableNames)
                            real_focal_vec(v) = T.real_focal_distance_mm(1);
                        else
                            warning('File %s missing real_focal_distance_mm.', filename);
                        end

                        if has_skull
                            if ismember('Isppa_skull', T.Properties.VariableNames)
                                maxIsppa_s_vec(v) = T.Isppa_skull(1);
                            else
                                warning('File %s missing Isppa_skull.', filename);
                            end

                            if ismember('riseT_skull', T.Properties.VariableNames)
                                riseT_s_vec(v) = T.riseT_skull(1);
                            else
                                warning('File %s missing riseT_skull.', filename);
                            end
                        end
                    else
                        warning('File not found: %s', filename);
                    end
                end

                results.sweep(:, i_cond)    = sweep_values;
                results.Ipa_target(:, i_cond) = isppa_vec;
                results.maxT(:, i_cond)     = maxT_vec - 37;
                Isppa_brain(:, i_cond)      = maxIsppa_b_vec;
                riseT_brain(:, i_cond)      = riseT_b_vec;
                real_focal_mm(:, i_cond)    = real_focal_vec;
                if has_skull
                    Isppa_skull(:, i_cond)  = maxIsppa_s_vec;
                    riseT_skull(:, i_cond)  = riseT_s_vec;
                end
            end
        end

        %% Plots

        tag = [sub_str '_' variant];

        plot_metric_with_patches(results.maxT, xlabels, 'max. Temp rise (deg C)', ...
            [tag '_temperature'], [], scatter_levels, scatter_colors, pn);

        plot_metric_with_patches(results.Ipa_target, xlabels, 'target intensity (W/cm^2)', ...
            [tag '_intensity'], [], scatter_levels, scatter_colors, pn);

        plot_metric_with_patches(Isppa_brain, xlabels, 'max Isppa brain (W/cm^2)', ...
            [tag '_maxIsppa_brain'], [], scatter_levels, scatter_colors, pn);

        plot_metric_with_patches(riseT_brain, xlabels, 'rise T brain (deg C)', ...
            [tag '_riseT_brain'], [], scatter_levels, scatter_colors, pn);

        plot_metric_with_patches(real_focal_mm, xlabels, 'real focal distance (mm)', ...
            [tag '_real_focal_distance_mm'], [], scatter_levels, scatter_colors, pn);

        if has_skull
            plot_metric_with_patches(Isppa_skull, xlabels, 'max Isppa skull (W/cm^2)', ...
                [tag '_maxIsppa_skull'], [], scatter_levels, scatter_colors, pn);

            plot_metric_with_patches(riseT_skull, xlabels, 'rise T skull (deg C)', ...
                [tag '_riseT_skull'], [], scatter_levels, scatter_colors, pn);
        end

        %% Difference plots (high - low sweep level)

        n_params = size(results.maxT, 2);
        param_idx = 1:n_params;

        h = figure;
        set(h, 'Units', 'normalized', 'Position', [0.25 0.25 0.25 0.3]);
        hold on;
        bar(param_idx, results.maxT(2,:) - results.maxT(1,:), ...
            'FaceColor', [.7 .0 .0], 'EdgeColor', [1 1 1]);
        set(gca, 'XTick', param_idx, 'XTickLabel', strrep(xlabels, '_', ' '));
        xtickangle(45);
        yline(0, '-k', 'LineWidth', 1);
        ylabel('Difference in max. Temp rise (deg C)');
        set(gca, 'Layer', 'top');
        h.Color = 'white'; h.InvertHardcopy = 'off';
        saveas(h, fullfile(pn.figures, [tag '_diff_temperature.png']));
        close(h);

        h = figure;
        set(h, 'Units', 'normalized', 'Position', [0.25 0.25 0.25 0.3]);
        hold on;
        bar(param_idx, results.Ipa_target(2,:) - results.Ipa_target(1,:), ...
            'FaceColor', [.7 .0 .0], 'EdgeColor', [1 1 1]);
        set(gca, 'XTick', param_idx, 'XTickLabel', strrep(xlabels, '_', ' '));
        xtickangle(45);
        yline(0, '-k', 'LineWidth', 1);
        ylabel('Difference in target intensity (W/cm^2)');
        set(gca, 'Layer', 'top');
        h.Color = 'white'; h.InvertHardcopy = 'off';
        saveas(h, fullfile(pn.figures, [tag '_diff_intensity.png']));
        close(h);

        close all;
    end
end
