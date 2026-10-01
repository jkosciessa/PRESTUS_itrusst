% PRESTUS_PARAMVAR_03_COMPARE_2D_3D
%
% Compares 2D axisymmetric and 3D simulation outputs for the two ITRUSST
% benchmarks:
%   Benchmark 1 — water/brain phantom:  sub-001 (2D)  vs  sub-101 (3D)
%   Benchmark 2 — skull/brain phantom:  sub-002 (2D)  vs  sub-102 (3D)
%
% For each benchmark the script:
%   1. Loads the output CSV tables and prints key scalar metrics side-by-side.
%   2. Loads the intensity Nifti, extracts the focal-axis profile (along the
%      beam propagation axis through the focus x/y voxel), and superimposes them.
%
% The 2D representative run used is the first available paramvar condition
% with all OTHER tissue parameters held at their mid-range values (affix
% _b_c_1 for sub-001/002).  The 3D run is the single baseline (no affix).
%
% See also: PRESTUS_PARAMVAR_01_BENCHMARKS_SIM,
%           PRESTUS_PARAMVAR_01_BENCHMARKS_3D_SIM

clear all; close all; clc;

% -------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
addpath('/project/2425122.01/v05/tools/nifti_toolbox');
pn.figures = fullfile(rootpath, 'figures', 'debug', 'compare_2D_3D');
if ~isfolder(pn.figures); mkdir(pn.figures); end

dir_2d = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a');
dir_3d = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a_3D');

% Resolution (mm per voxel) — same for both 2D and 3D configs
voxel_mm = 0.5;

% -------------------------------------------------------------------------
% Benchmark definitions
% -------------------------------------------------------------------------
benchmarks = struct( ...
    'label',     {'Benchmark 1 (water/brain)', 'Benchmark 2 (skull/brain)'}, ...
    'sub_2d',    {001, 002}, ...
    'sub_3d',    {101, 102}, ...
    'affix_2d',  {'_b_c_1', '_b_c_1'});   % representative 2D mid-baseline run

% Scalar metrics to compare
metric_fields  = {'Isppa', 'Isppa_brain', 'Ipa_target', 'real_focal_distance_mm', 'riseT_brain', 'maxT', 'MI_brain'};
metric_labels  = {'ISPPA global max (W/cm²)', 'ISPPA brain (W/cm²)', 'Ipa target (W/cm²)', 'Focal distance (mm)', ...
                  'Temp rise brain (°C)', 'Max temp (°C)', 'MI brain'};

% -------------------------------------------------------------------------
% Loop over benchmarks
% -------------------------------------------------------------------------
for b = 1:numel(benchmarks)

    bm       = benchmarks(b);
    sub2d    = bm.sub_2d;
    sub3d    = bm.sub_3d;
    affix_2d = bm.affix_2d;

    fprintf('\n========================================================\n');
    fprintf('  %s\n', bm.label);
    fprintf('  2D: sub-%03d  |  3D: sub-%03d\n', sub2d, sub3d);
    fprintf('========================================================\n');

    % ------------------------------------------------------------------
    % 1. Scalar metrics from CSV
    % ------------------------------------------------------------------
    csv_2d = fullfile(dir_2d, sprintf('sub-%03d', sub2d), ...
        sprintf('sub-%03d_phantom_output_table%s.csv', sub2d, affix_2d));
    csv_3d = fullfile(dir_3d, sprintf('sub-%03d', sub3d), ...
        sprintf('sub-%03d_phantom.csv', sub3d));

    ok_2d = isfile(csv_2d);
    ok_3d = isfile(csv_3d);

    if ~ok_2d; warning('2D CSV not found: %s', csv_2d); end
    if ~ok_3d; warning('3D CSV not found: %s', csv_3d); end

    if ok_2d && ok_3d
        T2d = readtable(csv_2d);
        T3d = readtable(csv_3d);

        fprintf('\n  %-30s  %12s  %12s  %12s\n', 'Metric', '2D (axi)', '3D (full)', 'Diff (%)');
        fprintf('  %s\n', repmat('-', 1, 70));
        for m = 1:numel(metric_fields)
            fld = metric_fields{m};
            lbl = metric_labels{m};
            v2d = NaN; v3d = NaN;
            if ismember(fld, T2d.Properties.VariableNames); v2d = T2d.(fld)(1); end
            if ismember(fld, T3d.Properties.VariableNames); v3d = T3d.(fld)(1); end
            if v2d == 37 && contains(fld, 'maxT'); v2d = NaN; end   % uninitialised sentinel
            pct = 100 * (v3d - v2d) / v2d;
            fprintf('  %-30s  %12.4f  %12.4f  %+11.1f%%\n', lbl, v2d, v3d, pct);
        end
    end

    % ------------------------------------------------------------------
    % 2. Focal-axis intensity profiles from Nifti
    % ------------------------------------------------------------------
    nii_2d = fullfile(dir_2d, sprintf('sub-%03d', sub2d), 'nii', ...
        sprintf('sub-%03d_phantom_T1w%s_intensity.nii.gz', sub2d, affix_2d));
    nii_3d = fullfile(dir_3d, sprintf('sub-%03d', sub3d), 'nii', ...
        sprintf('sub-%03d_phantom_T1w_intensity.nii.gz', sub3d));

    ok_nii_2d = isfile(nii_2d);
    ok_nii_3d = isfile(nii_3d);

    if ~ok_nii_2d; warning('2D Nifti not found: %s', nii_2d); end
    if ~ok_nii_3d; warning('3D Nifti not found: %s', nii_3d); end

    if ok_nii_2d && ok_nii_3d

        img2d = load_nii(nii_2d);
        img3d = load_nii(nii_3d);

        vol2d = double(img2d.img);
        vol3d = double(img3d.img);

        % Use the CSV brain-focus voxel for both 2D and 3D.
        % For 2D axisymmetric the cylindrical expansion places the beam
        % axis at the centre voxel, so Ix/Iy_brain_vox == centre for both.
        if ok_2d && ok_3d
            ix2d = T2d.Ix_brain_vox(1);
            iy2d = T2d.Iy_brain_vox(1);
            ix3d = T3d.Ix_brain_vox(1);
            iy3d = T3d.Iy_brain_vox(1);
        else
            ix2d = round(size(vol2d, 1) / 2);
            iy2d = round(size(vol2d, 2) / 2);
            ix3d = round(size(vol3d, 1) / 2);
            iy3d = round(size(vol3d, 2) / 2);
        end

        % Extract axial (z, dim 3) profile through focal x/y position
        prof_2d = squeeze(vol2d(ix2d, iy2d, :));
        prof_3d = squeeze(vol3d(ix3d, iy3d, :));

        % Remove negative values
        prof_2d(prof_2d < 0) = 0;
        prof_3d(prof_3d < 0) = 0;

        % Shift z so that the transducer face (z=1 mm from config) is at 0
        trans_z_mm = 1;   % trans_pos z = 2 voxels × 0.5 mm/vox = 1 mm
        z_2d = (0:numel(prof_2d)-1) * voxel_mm - trans_z_mm;
        z_3d = (0:numel(prof_3d)-1) * voxel_mm - trans_z_mm;

        % ---- plot ----
        figure('Name', bm.label, 'Color', 'w', 'Units', 'centimeters', ...
               'Position', [2 2 14 10]);

        plot(z_2d, prof_2d, 'b-',  'LineWidth', 1.5, 'DisplayName', '2D axisymmetric'); hold on;
        plot(z_3d, prof_3d, 'r--', 'LineWidth', 1.5, 'DisplayName', '3D full');
        xlabel('Distance from transducer face (mm)');
        ylabel('Intensity (W/cm²)');
        title(bm.label);
        legend('Location', 'best');
        xlim([0, max([z_2d(end), z_3d(end)])]);
        grid on; box on;

        % Save figure
        figname = fullfile(pn.figures, sprintf('compare_2D_3D_%s', ...
            lower(strrep(strrep(bm.label, ' ', '_'), '/', ''))));
        saveas(gcf, [figname, '.png']);
        saveas(gcf, [figname, '.eps'], 'epsc');
        fprintf('\n  Figure saved: %s.png/.eps\n', figname);
    end

end

fprintf('\nDone.\n');
