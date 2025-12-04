%% plot GM only

% Load the three NIfTI files
filename1 = '/Volumes/2425122.01/PRESTUS_itrusst/data/tussim/itrusst_protocol1/sub-001/sub-001_final_isppa_orig_coord_b_c_1.nii.gz';
filename2 = '/Volumes/2425122.01/PRESTUS_itrusst/data/tussim/itrusst_protocol1/sub-001/sub-001_final_isppa_orig_coord_b_c_2.nii.gz';
filename3 = '/Volumes/2425122.01/PRESTUS_itrusst/data/tussim/itrusst_protocol1/sub-001/sub-001_final_isppa_orig_coord_b_k_1.nii.gz';

V1 = niftiread(filename1);
V2 = niftiread(filename2);
V3 = niftiread(filename3);

% Extract central 2D matrices (central axial slice - assuming 3rd dim is axial)
[~, nx, nz] = size(V1);
mid_x = round(nx/2);
central_slice1 = squeeze(V1(mid_x, :, :));  % nz x ny
central_slice2 = squeeze(V2(mid_x, :, :));
central_slice3 = squeeze(V3(mid_x, :, :));

% Plot half of the 2nd dimension (ny/2, central row of 1st dimension)
center_row = round(size(central_slice1, 1)/2);
half_profile1 = central_slice1(center_row, :);
half_profile2 = central_slice2(center_row, :);
half_profile3 = central_slice3(center_row, :);

figure('Position', [100 100 600 250]);
hold on;
plot([1:140]-20,half_profile1, 'r-', 'LineWidth', 2, 'DisplayName', 'Low sound speed (b_c_1)');
plot([1:140]-20,half_profile2, 'r:', 'LineWidth', 2, 'DisplayName', 'High sound speed (b_c_2)');
plot([1:140]-20,half_profile3, 'k', 'LineWidth', 2, 'DisplayName', 'Avg. sound speed (b_k_1)');
hold off;

xlabel('Position (mm from bowl)');
ylabel('Intensity');
legend('Location', 'best');
grid on;
xlim([0 100])


%% add skull

% Load the three NIfTI files
filename1 = '/Volumes/2425122.01/PRESTUS_itrusst/data/tussim/itrusst_protocol1/sub-002/sub-002_final_isppa_orig_coord_b_c_1.nii.gz';
filename2 = '/Volumes/2425122.01/PRESTUS_itrusst/data/tussim/itrusst_protocol1/sub-002/sub-002_final_isppa_orig_coord_b_c_2.nii.gz';
filename3 = '/Volumes/2425122.01/PRESTUS_itrusst/data/tussim/itrusst_protocol1/sub-002/sub-002_final_isppa_orig_coord_b_k_1.nii.gz';

V1 = niftiread(filename1);
V2 = niftiread(filename2);
V3 = niftiread(filename3);

% Extract central 2D matrices (central axial slice - assuming 3rd dim is axial)
[~, nx, nz] = size(V1);
mid_x = round(nx/2);
central_slice1 = squeeze(V1(mid_x, :, :));  % nz x ny
central_slice2 = squeeze(V2(mid_x, :, :));
central_slice3 = squeeze(V3(mid_x, :, :));

% Plot half of the 2nd dimension (ny/2, central row of 1st dimension)
center_row = round(size(central_slice1, 1)/2);
half_profile1 = central_slice1(center_row, :);
half_profile2 = central_slice2(center_row, :);
half_profile3 = central_slice3(center_row, :);

% figure('Position', [100 100 600 250]);
hold on;
plot([1:140]-20,half_profile1, 'r-', 'LineWidth', 2);
plot([1:140]-20,half_profile2, 'r:', 'LineWidth', 2);
plot([1:140]-20,half_profile3, 'k', 'LineWidth', 2);
hold off;

xlabel('Position (mm from bowl)');
ylabel('Intensity');
legend('off'); grid on;
xlim([30 90])
set(findall(gcf,'-property','FontSize'),'FontSize',14)