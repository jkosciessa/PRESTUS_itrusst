% Set up root path and directories
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'))
rootpath = pwd;

% setup_list = {'medium_setup1', ...
%                 'medium_kPlan', ...
%                 'medium_reference', ...
%                 'medium_PRESTUS', ...
%                 'medium_BabelBrain', ...
%                 'medium_setup2'};

setup_list = {'medium_setup1', ...
                'medium_kPlan', ...
                'medium_PRESTUS', ...
                'medium_BabelBrain', ...
                'medium_setup2'};

subject_list = {'sub-101', 'sub-102'};

pn.figures = fullfile(rootpath, 'figures', 'pipeline');
if ~isfolder(pn.figures); mkdir(pn.figures); end

% Preallocate result matrices
n_setups = numel(setup_list);
n_subjects = numel(subject_list);

maxT_mat = nan(n_subjects, n_setups);
isppa_mat = nan(n_subjects, n_setups);       % target (Ipa_target)
isppa_brain_mat = nan(n_subjects, n_setups); % brain max (Isppa_brain)

% Loop over setups and subjects
for s = 1:n_setups
    setup = setup_list{s};
    for sub = 1:n_subjects
        subject = subject_list{sub};
        data_dir = fullfile(rootpath, 'data', 'tussim', setup, subject);
        filename = fullfile(data_dir, [subject '_phantom.csv']);
        if isfile(filename)
            T = readtable(filename);
            if all(ismember({'Ipa_target','Isppa_brain','maxT'}, T.Properties.VariableNames))
                maxT_mat(sub, s) = T.maxT(1)-37;
                isppa_mat(sub, s) = T.Ipa_target(1);
                isppa_brain_mat(sub, s) = T.Isppa_brain(1);
            else
                warning('File %s lacks required columns.', filename);
            end
        else
            warning('File %s not found.', filename);
        end
    end
end

% summary table (rows: subject x setup)
[sub_idx, set_idx] = ndgrid(1:n_subjects, 1:n_setups);
writetable(table(subject_list(sub_idx(:))', setup_list(set_idx(:))', maxT_mat(:), isppa_mat(:), isppa_brain_mat(:), ...
    'VariableNames', {'subject', 'setup', 'riseT_max', 'Ipa_target', 'Isppa_brain'}), ...
    fullfile(pn.figures, 'pipeline_summary.csv'));

% maxT_mat / isppa_mat / isppa_brain_mat: rows = subjects, cols = setups
% setup_list = {'setup 1', 'kPlan', 'reference', 'PRESTUS', 'BabelBrain', 'setup 2'};
setup_list = {'setup 1', 'kPlan', 'PRESTUS', 'BabelBrain', 'setup 2'};
subject_list = {'sub-101', 'sub-102'};

% ------- Maximum Temperature -------
h = figure;
for i = 1:2
    subplot(1,2,i);
    bar(maxT_mat(i,:), 'FaceColor', [.7 .2 .2]);
    set(gca, 'XTick', 1:size(isppa_mat,2), 'XTickLabel', setup_list, 'FontSize', 14);
    xlabel('Default Setup', 'FontSize', 14);
    ylabel('Temperature rise (°C)', 'FontSize', 14);
    title(['benchmark ', num2str(i)], 'FontSize', 14);
    ylim([min(maxT_mat(:))-0.5, max(maxT_mat(:))+0.5]);
    grid on;
    xtickangle(30);
    ylim([0, 3.5])
end
figureName = ['pipeline_temperature'];
h.Color = 'white';
h.InvertHardcopy = 'off';
saveas(h, fullfile(pn.figures, figureName), 'epsc');
saveas(h, fullfile(pn.figures, figureName), 'png');

% ------- Target Intensity -------
h = figure;
for i = 1:2
    subplot(1,2,i);
    hb = bar([isppa_mat(i,:); isppa_brain_mat(i,:)]');
    hb(1).FaceColor = [.2 .5 .7]; hb(2).FaceColor = [.6 .6 .6];
    set(gca, 'XTick', 1:size(isppa_mat,2), 'XTickLabel', setup_list, 'FontSize', 14);
    xlabel('Default Setup', 'FontSize', 14);
    ylabel('Intensity (W/cm^2)', 'FontSize', 14);
    if i == 1; legend(hb, {'target', 'brain max'}, 'Location', 'north', 'Orientation', 'horizontal', 'Box', 'off'); end
    title(['benchmark ', num2str(i)], 'FontSize', 14);
    grid on;
    xtickangle(30);
    ylim([0 60])
end
figureName = ['pipeline_intensity'];
h.Color = 'white';
h.InvertHardcopy = 'off';
saveas(h, fullfile(pn.figures, figureName), 'epsc');
saveas(h, fullfile(pn.figures, figureName), 'png');

% ------- Target Intensity -------
h = figure;
for i = 1:2
    subplot(1,2,i);
    hb = bar([isppa_mat(i,:); isppa_brain_mat(i,:)]');
    hb(1).FaceColor = [.2 .5 .7]; hb(2).FaceColor = [.6 .6 .6];
    set(gca, 'XTick', 1:size(isppa_mat,2), 'XTickLabel', setup_list, 'FontSize', 14);
    xlabel('Default Setup', 'FontSize', 14);
    ylabel('Intensity (W/cm^2)', 'FontSize', 14);
    if i == 1; legend(hb, {'target', 'brain max'}, 'Location', 'north', 'Orientation', 'horizontal', 'Box', 'off'); end
    title(['benchmark ', num2str(i)], 'FontSize', 14);
    grid on;
    xtickangle(30);
    ylim([0 20])
end
figureName = ['pipeline_intensity_rescaled'];
h.Color = 'white';
h.InvertHardcopy = 'off';
saveas(h, fullfile(pn.figures, figureName), 'epsc');
saveas(h, fullfile(pn.figures, figureName), 'png');