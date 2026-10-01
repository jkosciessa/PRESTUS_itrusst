if ~isfolder('/Volumes/2425122.01/PRESTUS_itrusst/figures/debug/heating_maps'); mkdir('/Volumes/2425122.01/PRESTUS_itrusst/figures/debug/heating_maps'); end

% Load the NIfTI files using niftiread (MATLAB R2017b+ built-in)
mask_file = '/Volumes/2425122.01/PRESTUS_itrusst/data/tussim/itrusst_protocol1/sub-002/sub-002_final_medium_masks_orig_coord_s_a0_1.nii.gz';
img1_file = '/Volumes/2425122.01/PRESTUS_itrusst/data/tussim/itrusst_protocol1/sub-002/sub-002_final_heating_orig_coord_s_a0_1.nii.gz';
img2_file = '/Volumes/2425122.01/PRESTUS_itrusst/data/tussim/itrusst_protocol1/sub-002/sub-002_final_heating_orig_coord_s_a0_2.nii.gz';

mask = niftiread(mask_file);
img1 = niftiread(img1_file);
img2 = niftiread(img2_file);

diff_img = img2-img1;  % High attenuation - low attenuation

% Extract central axial slice (adjust slice index if data is 3D; common for PRESTUS sims)
if ndims(mask) == 3
    mid_slice = round(size(mask,1)/2);
    mask_slice = squeeze(mask(mid_slice,:,:));
    diff_slice = squeeze(diff_img(mid_slice,:,:));
    img1_slice = squeeze(img1(mid_slice,:,:));
    img2_slice = squeeze(img2(mid_slice,:,:));
else
    mask_slice = mask;
    diff_slice = diff_img;
end

% Compute symmetric limits for zero-centered colormap
cmax = max(abs(diff_slice(:)));

% Create custom diverging colormap: blue (neg) - white (zero) - red (pos)
n = 128;
black_to_white = [linspace(0,1,n)' linspace(0,1,n)' linspace(0,1,n)'];
white_to_red = [ones(n,1) linspace(1,0,n)' linspace(1,0,n)'];
custom_cmap = [black_to_white; white_to_red];

% Plot
figure('Position', [100 100 600 400]);
hold on;
imagesc(diff_slice, [-cmax cmax]);  % Symmetric explicit limits
colormap(custom_cmap);  % Apply custom white-zero map
cb = colorbar;  % or cb = colorbar('eastoutside');
cb.Label.String = 'Heating Difference (^{\circ}C)';
cb.Label.FontSize = 14;
cb.Label.Interpreter = 'tex';  % or 'latex'hold on;
contour(mask_slice, [2.5, 2.5], 'k-', 'LineWidth', 1, 'LineStyle', '-');  % Dashed black mask contour
hold off;
axis image equal tight;
title('High vs low skull attenuation', 'FontSize', 16);
% xlabel('X (pixels)'); ylabel('Y (pixels)');
set(gca, 'FontSize', 14);
print(gcf, '/Volumes/2425122.01/PRESTUS_itrusst/figures/debug/heating_maps/attenuation_heating_diff.png', '-dpng', '-r300');
print(gcf, '/Volumes/2425122.01/PRESTUS_itrusst/figures/debug/heating_maps/attenuation_heating_diff.eps', '-depsc');

% Plot heating at low attenuation
custom_cmap = [white_to_red];
figure('Position', [100 100 600 400]);
hold on;
imagesc(img1_slice, [37 39]);
colormap(custom_cmap);  % Apply custom white-zero map
cb = colorbar;  % or cb = colorbar('eastoutside');
cb.Label.String = 'Temperature (^{\circ}C)';
cb.Label.FontSize = 14;
cb.Label.Interpreter = 'tex';  % or 'latex'hold on;
contour(mask_slice, [2.5, 2.5], 'k-', 'LineWidth', 1, 'LineStyle', '-');  % Dashed black mask contour
hold off;
axis image equal tight;
title('low attenuation', 'FontSize', 16);
% xlabel('X (pixels)'); ylabel('Y (pixels)');
set(gca, 'FontSize', 14);
print(gcf, '/Volumes/2425122.01/PRESTUS_itrusst/figures/debug/heating_maps/attenuation_heating_low.png', '-dpng', '-r300');
print(gcf, '/Volumes/2425122.01/PRESTUS_itrusst/figures/debug/heating_maps/attenuation_heating_low.eps', '-depsc');
