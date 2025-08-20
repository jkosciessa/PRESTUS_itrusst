% Set root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

% Updated paths for your data and figures
pn.figures = fullfile(rootpath, 'figures');
pn.data = fullfile(rootpath, 'data', 'tussim','itrusst_protocol1');

% Subjects 301 to 312 correspond to thickness 1 mm to 12 mm
subjects = 301:312;
thicknesses = 1:12; % Thickness in mm corresponding to subjects

% Preallocate vectors for intensity
brain_intensity = nan(length(subjects), 1);

% Loop over subjects to load brain intensity from CSV files
for i = 1:length(subjects)
    sub_id = subjects(i);
    thickness = thicknesses(i);
    
    % Construct path to CSV file for this subject
    % Assuming pattern: .../sub-303/sub-303_phantom_output_table.csv
    csv_file = fullfile(pn.data, sprintf('sub-%03d', sub_id), sprintf('sub-%03d_phantom_output_table.csv', sub_id));
    
    if isfile(csv_file)
        T = readtable(csv_file);
        % Check if 'isppa_at_target' exists in table
        if any(strcmp(T.Properties.VariableNames, 'isppa_at_target'))
            target_intensity(i) = T.isppa_at_target(1); % assume first row, brain intensity
            brain_intensity(i) = T.max_Isppa_brain(1); % assume first row, brain intensity
        else
            warning('Column ''isppa_at_target'' not found in %s', csv_file);
        end
    else
        warning('File %s not found.', csv_file);
    end
end

%% Plot target intensity as a function of thickness

h = figure;
set(h, 'Units', 'normalized');
set(h, 'Position', [0.25 0.25 0.15 0.2]);
p1 = plot(thicknesses, target_intensity, 'k-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
p2 = plot(thicknesses, brain_intensity, '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
grid on;
% Add vertical red line at x = 6.5
xline(6.5, 'r', 'LineWidth', 2);
xlabel('Skull Thickness (mm)');
ylabel('Intensity (W/cm2)'); % Adjust label according to intensity measure
xlim([1, 12]);
legend([p1, p2], {'target'; 'brain'})
legend('boxoff')
set(gca, 'FontSize', 14);

% Save plot
if ~exist(pn.figures, 'dir')
    mkdir(pn.figures);
end
figureName = 'target_intensity_vs_thickness';
saveas(h, fullfile(pn.figures, figureName), 'png');
saveas(h, fullfile(pn.figures, figureName), 'epsc');