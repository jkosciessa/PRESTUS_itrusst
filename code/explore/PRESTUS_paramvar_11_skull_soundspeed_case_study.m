% PRESTUS_PARAMVAR_11_SKULL_SOUNDSPEED_CASE_STUDY
%
% Case-study figures for skull sound_speed (c_{s}), the OAT parameter with
% the largest impact on real focal distance and -3dB focal volume in the
% PRESTUS_paramvar_06_impact_direction_plot.m tornado plots — Benchmark 2
% (sub-002, skull + brain).
%
% (a) Axial acoustic-intensity and thermal-heatrise profiles for the bottom
%     (low c_skull) vs top (high c_skull) sweep bounds, superimposed, with
%     the peak location marked and the CSV-reported focal distance / -3dB
%     focal volume annotated for cross-reference against the tornado plots.
%
% (b) A 2D difference map of the -3dB focal-volume mask (focal_mask_3dB,
%     already saved per condition by the PRESTUS toolbox — see
%     compute_focal_metrics.m), showing WHERE the focal zone is gained,
%     lost, or shared as skull sound speed moves from its bottom to its top
%     literature bound. Uses a 4-category map (neither / low-c only /
%     high-c only / both) rather than a raw difference, since a plain
%     high-minus-low subtraction cannot distinguish "outside the focal zone
%     in both conditions" from "inside the focal zone in both conditions"
%     (both give zero).
%
% Note: the 2D maps are the *_intensity.nii.gz / *_focal_mask_3dB.nii.gz
% NIfTIs already exported per condition (a 2-D axisymmetric cross-section
% cropped to the active simulation domain). The CSV's focal_vol_3dB_mm3 is
% computed on the full internal simulation grid (including the PML margin
% this export omits), so the two will not numerically match voxel-for-voxel
% — the map is a spatial illustration of the reported CSV totals, not a
% pixel-accurate reproduction of them.
%
% Data:   data/tussim/itrusst_protocol1a/sub-002 (2D axisymmetric, 0.5 mm grid)
% Output: figures/explore/skull_soundspeed/
%
% See also: PRESTUS_PARAMVAR_04_AXIAL_PROFILES_2D,
%           PRESTUS_PARAMVAR_06_IMPACT_DIRECTION_PLOT, COMPUTE_FOCAL_METRICS

clear; close all; clc;

%% ── Paths ───────────────────────────────────────────────────────────────────

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));

pn.data    = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a', 'sub-002');
pn.nii     = fullfile(pn.data, 'nii');
pn.figures = fullfile(rootpath, 'figures', 'explore', 'skull_soundspeed');
if ~isfolder(pn.figures); mkdir(pn.figures); end

%% ── Configuration: skull sound_speed OAT sweep (Benchmark 2) ────────────────

tissue_suffix = '_s';    % skull
prop_suffix   = '_c_';   % sound_speed
levels        = [1, 2];  % 1 = bottom (low c_skull), 2 = top (high c_skull)
level_labels  = {'low c_{skull} (bottom bound)', 'high c_{skull} (top bound)'};
level_colors  = {[0.2 0.5 0.9], [0.9 0.3 0.2]};   % blue = low, red = high

voxel_mm  = 0.5;
trans_vox = 1;   % transducer face at z-voxel 1 (config_itrusst_protocol.yaml)

affix = @(lv) sprintf('%s%s%d', tissue_suffix, prop_suffix, lv);

%% ── Load per-level CSV summary metrics (for titles / cross-reference) ───────

csv_metrics = struct('real_focal_distance_mm', {}, 'focal_vol_3dB_mm3', {}, 'Isppa_brain', {});
for i = 1:numel(levels)
    csv_file = fullfile(pn.data, sprintf('sub-002_phantom%s.csv', affix(levels(i))));
    if ~isfile(csv_file)
        error('CSV not found: %s', csv_file);
    end
    T = readtable(csv_file);
    csv_metrics(i).real_focal_distance_mm = T.real_focal_distance_mm(1);
    csv_metrics(i).focal_vol_3dB_mm3      = T.focal_vol_3dB_mm3(1);
    csv_metrics(i).Isppa_brain            = T.Isppa_brain(1);
end

fprintf('Skull sound_speed OAT sweep — sub-002 (skull + brain)\n');
fprintf('  %-28s  %10s  %10s\n', 'metric', level_labels{1}, level_labels{2});
fprintf('  %-28s  %10.2f  %10.2f\n', 'real_focal_distance_mm', csv_metrics(1).real_focal_distance_mm, csv_metrics(2).real_focal_distance_mm);
fprintf('  %-28s  %10.2f  %10.2f\n', 'focal_vol_3dB_mm3',      csv_metrics(1).focal_vol_3dB_mm3,      csv_metrics(2).focal_vol_3dB_mm3);
fprintf('  %-28s  %10.2f  %10.2f\n', 'Isppa_brain',            csv_metrics(1).Isppa_brain,            csv_metrics(2).Isppa_brain);

%% ── (a) Axial profile comparison: acoustic + thermal, low vs high c_skull ──

modalities = struct( ...
    'name',   {'acoustic',            'thermal'}, ...
    'suffix', {'_intensity.nii.gz',   '_heatrise.nii.gz'}, ...
    'ylabel', {'Intensity (W/cm^2)',  '\DeltaT (^{\circ}C)'});

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [1 1 26 11], ...
             'Name', 'sub-002 — skull sound_speed — axial profiles');

for mo = 1:numel(modalities)
    mod = modalities(mo);
    ax  = subplot(1, numel(modalities), mo);
    hold(ax, 'on');

    for i = 1:numel(levels)
        nii_path = fullfile(pn.nii, sprintf('sub-002_phantom_T1w%s%s', affix(levels(i)), mod.suffix));
        if ~isfile(nii_path)
            warning('File not found: %s', nii_path);
            continue;
        end

        vol  = double(niftiread(nii_path));
        cx   = round(size(vol, 1) / 2);
        nz   = size(vol, 2);
        z_mm = ((1:nz) - trans_vox) * voxel_mm;
        prof = vol(cx, :);
        prof(prof < 0) = 0;

        plot(ax, z_mm, prof, 'Color', level_colors{i}, 'LineWidth', 2, ...
             'DisplayName', level_labels{i});

        [~, pk_idx] = max(prof);
        xline(ax, z_mm(pk_idx), '--', 'Color', level_colors{i}, 'LineWidth', 1, ...
              'HandleVisibility', 'off');
    end

    xlabel(ax, 'Distance from transducer (mm)');
    ylabel(ax, mod.ylabel, 'Interpreter', 'tex');
    title(ax, sprintf('%s axial profile', mod.name));
    grid(ax, 'on'); box(ax, 'on');
    if mo == 1
        legend(ax, 'Location', 'best');
    end
end

sgtitle(sprintf(['sub-002 (skull+brain), c_{skull}: bottom vs top literature bound\n' ...
    'real focal distance: %.1f \\rightarrow %.1f mm   |   -3dB focal volume: %.0f \\rightarrow %.0f mm^3'], ...
    csv_metrics(1).real_focal_distance_mm, csv_metrics(2).real_focal_distance_mm, ...
    csv_metrics(1).focal_vol_3dB_mm3, csv_metrics(2).focal_vol_3dB_mm3), 'Interpreter', 'tex');

figname = fullfile(pn.figures, 'sub-002_skull_soundspeed_axial_profiles');
saveas(fig, [figname '.png']);
saveas(fig, [figname '.eps'], 'epsc');
fprintf('\nSaved: %s.png/.eps\n', figname);

%% ── (b) 2D map of focal-volume (-3dB mask) change: low vs high c_skull ─────

mask_files = cell(1, numel(levels));
for i = 1:numel(levels)
    mask_files{i} = fullfile(pn.nii, sprintf('sub-002_phantom_T1w%s_focal_mask_3dB.nii.gz', affix(levels(i))));
    if ~isfile(mask_files{i})
        error('Focal mask NIfTI not found: %s', mask_files{i});
    end
end

mask_lo = double(niftiread(mask_files{1})) > 0;   % low c_skull
mask_hi = double(niftiread(mask_files{2})) > 0;   % high c_skull

% 4-category map: 0 = neither, 1 = low-c only, 2 = high-c only, 3 = both.
% A raw high-minus-low difference cannot tell "outside the focal zone in
% both conditions" apart from "inside the focal zone in both conditions"
% (both give zero), so use categories instead.
category = zeros(size(mask_lo));
category(mask_lo & ~mask_hi) = 1;   % lost as c_skull increases
category(~mask_lo & mask_hi) = 2;   % gained as c_skull increases
category(mask_lo & mask_hi)  = 3;   % retained in both

nx   = size(category, 1);
nz   = size(category, 2);
z_mm = ((1:nz) - trans_vox) * voxel_mm;
x_mm = ((1:nx) - round(nx / 2)) * voxel_mm;   % centred on the beam axis

cat_colors = [ ...
    0.93 0.93 0.93;   % 0 neither
    0.20 0.50 0.90;   % 1 low-c only  (lost at high c)
    0.90 0.30 0.20;   % 2 high-c only (gained at high c)
    0.55 0.35 0.75];  % 3 both

fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [1 1 18 14], ...
              'Name', 'sub-002 — skull sound_speed — focal volume change');
imagesc(z_mm, x_mm, category, [0 3]);
colormap(gca, cat_colors);
set(gca, 'YDir', 'normal');
axis image;

cb = colorbar('Ticks', 0:3, 'TickLabels', ...
    {'neither', 'low c only (lost)', 'high c only (gained)', 'both (retained)'});
cb.Label.String = '-3dB focal zone membership';

xlabel('Distance from transducer (mm)');
ylabel('Lateral distance from beam axis (mm)');
title(sprintf(['sub-002 (skull+brain), c_{skull}: -3dB focal-zone change, low \\rightarrow high bound\n' ...
    'reported -3dB focal volume: %.0f \\rightarrow %.0f mm^3 (full simulation grid, CSV)'], ...
    csv_metrics(1).focal_vol_3dB_mm3, csv_metrics(2).focal_vol_3dB_mm3), 'Interpreter', 'tex');
set(gca, 'FontSize', 11, 'Layer', 'top');

figname2 = fullfile(pn.figures, 'sub-002_skull_soundspeed_focal_volume_change');
saveas(fig2, [figname2 '.png']);
saveas(fig2, [figname2 '.eps'], 'epsc');
fprintf('Saved: %s.png/.eps\n', figname2);

fprintf('\nDone.\n');
