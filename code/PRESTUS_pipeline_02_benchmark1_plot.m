% Set up root path and directories
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'))
rootpath = pwd;

setup_list = {'itrusst_protocol1_bottom', 'itrusst_protocol1_mid', ...
              'itrusst_protocol1_PRESTUS', 'itrusst_protocol1_kPlan', ...
              'itrusst_protocol1_BabelBrain', 'itrusst_protocol1_top'};

subject_list = {'sub-001', 'sub-002'};

pn.figures = fullfile(rootpath, 'figures');

% Preallocate result matrices
n_setups = numel(setup_list);
n_subjects = numel(subject_list);

maxT_mat = nan(n_subjects, n_setups);
isppa_mat = nan(n_subjects, n_setups);

% Loop over setups and subjects
for s = 1:n_setups
    setup = setup_list{s};
    for sub = 1:n_subjects
        subject = subject_list{sub};
        data_dir = fullfile(rootpath, 'data', 'tussim', setup, subject);
        filename = fullfile(data_dir, [subject '_phantom_output_table.csv']);
        if isfile(filename)
            T = readtable(filename);
            if all(ismember({'isppa_at_target','maxT'}, T.Properties.VariableNames))
                maxT_mat(sub, s) = T.maxT(1);
                isppa_mat(sub, s) = T.isppa_at_target(1);
            else
                warning('File %s lacks required columns.', filename);
            end
        else
            warning('File %s not found.', filename);
        end
    end
end

% Assume maxT_mat and isppa_mat are [2 x 4] arrays: rows=subjects, cols=setups
setup_list = {'bottom', 'mid', 'PRESTUS', 'kPlan', 'BabelBrain', 'top'};
subject_list = {'sub-001', 'sub-002'};

% ------- Maximum Temperature -------
h = figure;
for i = 1:2
    subplot(1,2,i);
    bar(maxT_mat(i,:), 'FaceColor', [.2 .5 .7]);
    set(gca, 'XTick', 1:size(isppa_mat,2), 'XTickLabel', setup_list, 'FontSize', 14);
    xlabel('Default Setup', 'FontSize', 14);
    ylabel('Maximum Temperature (°C)', 'FontSize', 14);
    title(['benchmark ', num2str(i)], 'FontSize', 14);
    ylim([min(maxT_mat(:))-0.5, max(maxT_mat(:))+0.5]);
    grid on;
    xtickangle(30);
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
    bar(isppa_mat(i,:), 'FaceColor', [.7 .2 .2]);
    set(gca, 'XTick', 1:size(isppa_mat,2), 'XTickLabel', setup_list, 'FontSize', 14);
    xlabel('Default Setup', 'FontSize', 14);
    ylabel('Target Intensity', 'FontSize', 14);
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