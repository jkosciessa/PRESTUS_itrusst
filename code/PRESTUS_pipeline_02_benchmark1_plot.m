% Pipeline comparison (PRESTUS_pipeline_01_benchmark1_sim.m): bar plots per
% benchmark of
%   - heating at the focus: temperature rise at the target voxel
%     (focus_pos_vox) of the heat-rise NIfTI (beam-axis plane); heating
%     setups Temp-H / Temp-L and tools
%   - maximum brain heating, riseT_brain: Temp-H / Temp-L and tools; with a
%     skull this maximum lies next to the skull, not at the focus
%   - skull heating, riseT_skull (skull benchmark only): Temp-H / Temp-L and tools
%   - maximum brain intensity (Isppa_brain): transmission setups Trans-H /
%     Trans-L (acoustic-only) and tools
% Output: figures/pipeline/pipeline_{temperature_focus,temperature_brain,temperature_skull,
%         intensity,intensity_rescaled}.{eps,png}, pipeline_summary.csv

% Set up root path and directories
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'))
rootpath = pwd;

tool_dirs   = {'medium_kPlan', 'medium_PRESTUS', 'medium_BabelBrain'};
tool_labels = {'kPlan', 'PRESTUS', 'BabelBrain'};
heat_dirs   = [{'medium_heat_max'}, tool_dirs, {'medium_heat_min'}];
heat_labels = [{'Temp-H'}, tool_labels, {'Temp-L'}];
trans_dirs   = [{'medium_trans_max'}, tool_dirs, {'medium_trans_min'}];
trans_labels = [{'Trans-H'}, tool_labels, {'Trans-L'}];

subject_list = {'sub-101', 'sub-102'};
benchmark_titles = {'no skull (sub-101)', 'skull (sub-102)'};

pn.figures = fullfile(rootpath, 'figures', 'pipeline');
if ~isfolder(pn.figures); mkdir(pn.figures); end

% rows = subjects, cols = setups
riseT_focus = read_focus_rise(rootpath, heat_dirs, subject_list);
riseT_brain = read_metric(rootpath, heat_dirs, subject_list, 'riseT_brain');
riseT_skull = read_metric(rootpath, heat_dirs, subject_list, 'riseT_skull');
ipa_target  = read_metric(rootpath, trans_dirs, subject_list, 'Ipa_target');
isppa_brain = read_metric(rootpath, trans_dirs, subject_list, 'Isppa_brain');

% summary table (rows: subject x setup; thermal only for the heating
% setups and tools, intensity only for the transmission setups and tools)
setup_dirs = unique([heat_dirs, trans_dirs], 'stable');
summary = table();
for sub = 1:numel(subject_list)
    for s = 1:numel(setup_dirs)
        i_heat = find(strcmp(heat_dirs, setup_dirs{s}));
        i_trans = find(strcmp(trans_dirs, setup_dirs{s}));
        row = table(subject_list(sub), setup_dirs(s), NaN, NaN, NaN, NaN, NaN, ...
            'VariableNames', {'subject', 'setup', 'riseT_focus', 'riseT_brain', 'riseT_skull', 'Ipa_target', 'Isppa_brain'});
        if ~isempty(i_heat)
            row.riseT_focus = riseT_focus(sub, i_heat);
            row.riseT_brain = riseT_brain(sub, i_heat);
            row.riseT_skull = riseT_skull(sub, i_heat);
        end
        if ~isempty(i_trans)
            row.Ipa_target = ipa_target(sub, i_trans);
            row.Isppa_brain = isppa_brain(sub, i_trans);
        end
        summary = [summary; row]; %#ok<AGROW>
    end
end
writetable(summary, fullfile(pn.figures, 'pipeline_summary.csv'));

% ------- Temperature rise at the focus -------
h = figure;
for i = 1:2
    subplot(1,2,i);
    bar(riseT_focus(i,:), 'FaceColor', [.7 .2 .2]);
    set(gca, 'XTick', 1:numel(heat_labels), 'XTickLabel', heat_labels, 'FontSize', 14);
    xlabel('Setup', 'FontSize', 14);
    ylabel('Temperature rise at focus (°C)', 'FontSize', 14);
    title(benchmark_titles{i}, 'FontSize', 14);
    grid on;
    xtickangle(30);
    ylim([0, ceil(2 * max(riseT_focus(:), [], 'omitnan')) / 2]);
end
save_figure(h, pn.figures, 'pipeline_temperature_focus');

% ------- Maximum brain temperature rise -------
h = figure;
for i = 1:2
    subplot(1,2,i);
    bar(riseT_brain(i,:), 'FaceColor', [.7 .2 .2]);
    set(gca, 'XTick', 1:numel(heat_labels), 'XTickLabel', heat_labels, 'FontSize', 14);
    xlabel('Setup', 'FontSize', 14);
    ylabel('Brain temperature rise (°C)', 'FontSize', 14);
    title(benchmark_titles{i}, 'FontSize', 14);
    grid on;
    xtickangle(30);
    ylim([0, ceil(2 * max(riseT_brain(:), [], 'omitnan')) / 2]);
end
save_figure(h, pn.figures, 'pipeline_temperature_brain');

% ------- Skull temperature rise (skull benchmark only) -------
h = figure;
i = find(strcmp(subject_list, 'sub-102'));
bar(riseT_skull(i,:), 'FaceColor', [.7 .2 .2]);
set(gca, 'XTick', 1:numel(heat_labels), 'XTickLabel', heat_labels, 'FontSize', 14);
xlabel('Setup', 'FontSize', 14);
ylabel('Skull temperature rise (°C)', 'FontSize', 14);
title(benchmark_titles{i}, 'FontSize', 14);
grid on;
xtickangle(30);
ylim([0, ceil(2 * max(riseT_skull(i,:), [], 'omitnan')) / 2]);
save_figure(h, pn.figures, 'pipeline_temperature_skull');

% ------- Target Intensity (full and rescaled y-axis) -------
names = {'pipeline_intensity', 'pipeline_intensity_rescaled'};
ylimits = {[0 60], [0 30]};
for k = 1:2
    h = figure;
    for i = 1:2
        subplot(1,2,i);
        bar(isppa_brain(i,:), 'FaceColor', [.6 .6 .6]);
        set(gca, 'XTick', 1:numel(trans_labels), 'XTickLabel', trans_labels, 'FontSize', 14);
        xlabel('Setup', 'FontSize', 14);
        ylabel('Brain I_{SPPA} (W/cm^2)', 'FontSize', 14);
        title(benchmark_titles{i}, 'FontSize', 14);
        grid on;
        xtickangle(30);
        ylim(ylimits{k})
    end
    save_figure(h, pn.figures, names{k});
end

%% local functions

function M = read_metric(rootpath, setup_dirs, subject_list, metric)
% metric per subject (rows) x setup (cols) from <setup>/<subject>/<subject>_phantom.csv
M = nan(numel(subject_list), numel(setup_dirs));
for s = 1:numel(setup_dirs)
    for sub = 1:numel(subject_list)
        subject = subject_list{sub};
        filename = fullfile(rootpath, 'data', 'tussim', setup_dirs{s}, subject, [subject '_phantom.csv']);
        if ~isfile(filename)
            warning('File %s not found.', filename);
            continue
        end
        T = readtable(filename);
        if ismember(metric, T.Properties.VariableNames)
            M(sub, s) = T.(metric)(1);
        else
            warning('File %s lacks %s.', filename, metric);
        end
    end
end
end

function M = read_focus_rise(rootpath, setup_dirs, subject_list)
% temperature rise at the target voxel per subject (rows) x setup (cols):
% <subject>_phantom_T1w_heatrise.nii.gz holds the beam-axis plane
% [lateral x axial]; the target is (focus_pos_vox_1, focus_pos_vox_3) of the CSV
M = nan(numel(subject_list), numel(setup_dirs));
for s = 1:numel(setup_dirs)
    for sub = 1:numel(subject_list)
        subject = subject_list{sub};
        sim_dir = fullfile(rootpath, 'data', 'tussim', setup_dirs{s}, subject);
        csv_file = fullfile(sim_dir, [subject '_phantom.csv']);
        nii_file = fullfile(sim_dir, 'nii', [subject '_phantom_T1w_heatrise.nii.gz']);
        if ~isfile(csv_file) || ~isfile(nii_file)
            warning('Heat-rise map or CSV missing in %s.', sim_dir);
            continue
        end
        T = readtable(csv_file);
        rise = double(niftiread(nii_file));
        M(sub, s) = rise(T.focus_pos_vox_1(1), T.focus_pos_vox_3(1));
    end
end
end

function save_figure(h, out_dir, name)
h.Color = 'white';
h.InvertHardcopy = 'off';
saveas(h, fullfile(out_dir, name), 'epsc');
saveas(h, fullfile(out_dir, name), 'png');
end
