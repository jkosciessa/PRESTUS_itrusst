% Set root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

% Updated paths for your data and figures
pn.figures = fullfile(rootpath, 'figures', 'thickness');
if ~isfolder(pn.figures); mkdir(pn.figures); end
pn.data = fullfile(rootpath, 'data', 'tussim','itrusst_protocol');

% sub-301..312 = 1..12 mm, sub-351..361 = 1.5..11.5 mm
addpath(fullfile(rootpath, 'code', 'functions'));
[subjects, thicknesses] = itrusst_thickness_subjects(1:0.5:12);

% Preallocate vectors for intensity
brain_intensity = nan(length(subjects), 1);
target_intensity = nan(length(subjects), 1);

% Loop over subjects to load brain intensity from CSV files
for i = 1:length(subjects)
    sub_id = subjects(i);
    thickness = thicknesses(i);
    
    % Construct path to CSV file for this subject
    % Assuming pattern: .../sub-303/sub-303_phantom_output_table.csv
    csv_file = fullfile(pn.data, sprintf('sub-%03d', sub_id), sprintf('sub-%03d_phantom.csv', sub_id));
    
    if isfile(csv_file)
        T = readtable(csv_file);
        % Check if 'Ipa_target' exists in table
        if any(strcmp(T.Properties.VariableNames, 'Ipa_target'))
            target_intensity(i) = T.Ipa_target(1); % assume first row, brain intensity
            brain_intensity(i) = T.Isppa_brain(1); % assume first row, brain intensity
        else
            warning('Column ''Ipa_target'' not found in %s', csv_file);
        end
    else
        warning('File %s not found.', csv_file);
    end
end

%% Plot target and brain intensity as a function of thickness

h = figure;
set(gcf,'renderer','Painters');
set(h, 'Units', 'normalized');
set(h, 'Position', [0.25 0.25 0.15 0.2]);
p1 = plot(thicknesses, target_intensity, 'k-o', 'LineWidth', 2, 'MarkerSize', 4);
hold on;
p2 = plot(thicknesses, brain_intensity, '-o', 'LineWidth', 2, 'MarkerSize', 4, 'Color', [0.5 0.5 0.5]);
grid on;
% Add vertical red line at x = 6.5
xline(6.5, 'r', 'LineWidth', 2);
% Add vertical blue line at x = 3
xline(3, 'b', 'LineWidth', 2);
% Add vertical mean intensity
% yline(mean(target_intensity), 'r', 'LineWidth', 1, 'LineStyle', ':');
xlabel('Skull Thickness (mm)');
ylabel('Intensity (W/cm^2)');
xlim([1, 12]);
legend([p1, p2], {'target (Ipa\_target)'; 'brain max (Isppa\_brain)'}, 'Location', 'best')
legend('boxoff')
set(gca, 'FontSize', 12);

% Save plot
figureName = 'intensity_vs_thickness';
saveas(h, fullfile(pn.figures, figureName), 'png');
saveas(h, fullfile(pn.figures, figureName), 'epsc');
% Save table
writetable(table(subjects(:), thicknesses(:), target_intensity, brain_intensity, ...
    'VariableNames', {'subject_id', 'thickness_mm', 'Ipa_target', 'Isppa_brain'}), ...
    fullfile(pn.figures, 'intensity_vs_thickness.csv'));
