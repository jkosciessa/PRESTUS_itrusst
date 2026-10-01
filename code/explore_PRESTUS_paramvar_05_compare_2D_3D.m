% PRESTUS_PARAMVAR_05_COMPARE_2D_3D
%
% Benchmark: compares 2D axisymmetric and 3D simulation outputs for the
% high-resolution (0.25 mm) ITRUSST phantoms:
%
%   Benchmark 1 — water/brain:  sub-201 (2D axi)  vs  sub-201 (3D full)
%   Benchmark 2 — skull/brain:  sub-202 (2D axi)  vs  sub-202 (3D full)
%
% The "representative" 2D run is the mid-baseline condition (_b_c_1):
% brain sound speed at its bottom value, all other parameters at mid.
%
% For each benchmark the script:
%   1. Extracts the focal-axis (z) intensity profile from the center voxel
%      of each simulation volume.
%   2. Extracts the focal-axis (z) heatrise profile when available
%      (_heatrise.nii.gz — requires thermal sims to have completed).
%   3. Superimposes 2D and 3D profiles (absolute and peak-normalised) for
%      each modality.
%   4. Quantifies correspondence:
%        - Pearson correlation of the two profiles
%        - Peak ratio (3D / 2D)
%        - Focal distance difference (mm)
%        - FWHM difference (mm)
%
% 2D data:  data/tussim/itrusst_protocol1a/sub-201|202  (existing)
% 3D data:  data/tussim/itrusst_protocol1a_3D/sub-201|202  (run when available)
%
% If the 3D data directory is absent the script errors with a clear message.
% Thermal panels are skipped with a warning if _heatrise.nii.gz is missing
% (thermal sims may still be running — check log for 400/400 pulse trains).
%
% See also: PRESTUS_PARAMVAR_04_AXIAL_PROFILES_2D,
%           PRESTUS_PARAMVAR_03_COMPARE_2D_3D

clear; close all; clc;

% -------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code'));
addpath('/project/2425122.01/v05/tools/nifti_toolbox');
if isfolder('/Volumes/2425122.01/v05/tools/nifti_toolbox')
    addpath('/Volumes/2425122.01/v05/tools/nifti_toolbox');
end

pn.figures = fullfile(rootpath, 'figures', 'itrusst_protocol1a', 'compare_2D_3D');
if ~isfolder(pn.figures); mkdir(pn.figures); end

dir_2d = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a');
dir_3d = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a_3D');

if ~isfolder(dir_3d)
    error(['3D output directory not found:\n  %s\n' ...
           'Run PRESTUS_PARAMVAR_01_BENCHMARKS_3D_SIM for sub-201/202 first.'], dir_3d);
end

% -------------------------------------------------------------------------
% Grid parameters (must match config_itrusst_protocol1_025mm.yaml)
% -------------------------------------------------------------------------
voxel_mm_2d  = 0.25;
trans_vox_2d = 4;     % z-voxel of transducer face (2D axi config: trans_pos z = 4)

voxel_mm_3d  = 0.25;
trans_vox_3d = 4;     % adjust if 3D config differs

% Representative 2D run (brain c at bottom, all else mid)
affix_2d_rep = '_b_c_1';

% -------------------------------------------------------------------------
% Benchmark definitions
% -------------------------------------------------------------------------
benchmarks = struct( ...
    'label',   {'Benchmark 1 (water/brain)', 'Benchmark 2 (skull/brain)'}, ...
    'sub_id',  {201, 202});

% -------------------------------------------------------------------------
% Modality definitions
% -------------------------------------------------------------------------
modalities = struct( ...
    'name',   {'acoustic',          'thermal'}, ...
    'suffix', {'_intensity.nii.gz', '_heatrise.nii.gz'}, ...
    'ylabel', {'Intensity (W/cm²)', '\DeltaT (°C)'});

% -------------------------------------------------------------------------
% Loop over benchmarks
% -------------------------------------------------------------------------
results_table = table();

for b = 1:numel(benchmarks)

    bm      = benchmarks(b);
    sub_id  = bm.sub_id;
    sub_str = sprintf('sub-%03d', sub_id);

    fprintf('\n============================================================\n');
    fprintf('  %s  (%s)\n', sub_str, bm.label);
    fprintf('============================================================\n');

    nii_2d_dir = fullfile(dir_2d, sub_str, 'nii');
    nii_3d_dir = fullfile(dir_3d, sub_str, 'nii');

    row_data = {sub_id, bm.label};
    row_vars = {'sub_id', 'label'};

    for mo = 1:numel(modalities)
        mod = modalities(mo);

        nii_2d_path = fullfile(nii_2d_dir, ...
            sprintf('%s_phantom_T1w%s%s', sub_str, affix_2d_rep, mod.suffix));
        nii_3d_path = fullfile(nii_3d_dir, ...
            sprintf('%s_phantom_T1w%s', sub_str, mod.suffix));

        ok_2d = isfile(nii_2d_path);
        ok_3d = isfile(nii_3d_path);

        if ~ok_2d || ~ok_3d
            if strcmp(mod.name, 'acoustic')
                if ~ok_2d; warning('2D acoustic NII not found: %s', nii_2d_path); end
                if ~ok_3d; warning('3D acoustic NII not found: %s\nRun 3D sim first.', nii_3d_path); end
                break;
            else
                fprintf('  [thermal] NII not yet available — skipping.\n');
                fprintf('  (thermal sims may be running; check log for 400/400 pulse trains)\n');
                continue;
            end
        end

        % Load and extract profiles
        img2d = load_nii(nii_2d_path);  vol2d = double(img2d.img);  vol2d(vol2d < 0) = 0;
        img3d = load_nii(nii_3d_path);  vol3d = double(img3d.img);  vol3d(vol3d < 0) = 0;

        cx2d = round(size(vol2d, 1) / 2);  cy2d = round(size(vol2d, 2) / 2);
        cx3d = round(size(vol3d, 1) / 2);  cy3d = round(size(vol3d, 2) / 2);

        z_2d = ((1:size(vol2d, 3)) - trans_vox_2d) * voxel_mm_2d;
        z_3d = ((1:size(vol3d, 3)) - trans_vox_3d) * voxel_mm_3d;

        prof2d = squeeze(vol2d(cx2d, cy2d, :));
        prof3d = squeeze(vol3d(cx3d, cy3d, :));

        m2d = profile_metrics(prof2d, z_2d);
        m3d = profile_metrics(prof3d, z_3d);

        % Interpolate to common z-axis for correlation
        z_com = max(z_2d(1), z_3d(1)) : voxel_mm_2d : min(z_2d(end), z_3d(end));
        rmat  = corrcoef(interp1(z_2d, prof2d, z_com, 'linear', 0), ...
                         interp1(z_3d, prof3d, z_com, 'linear', 0));
        r     = rmat(1, 2);

        peak_ratio    = m3d.peak_val    / m2d.peak_val;
        focal_diff_mm = m3d.focal_dist_mm - m2d.focal_dist_mm;
        fwhm_diff_mm  = m3d.fwhm_mm      - m2d.fwhm_mm;

        fprintf('\n  [%s]\n', upper(mod.name));
        fprintf('  %-28s  %10s  %10s  %10s\n', 'Metric', '2D axi', '3D full', 'Diff / ratio');
        fprintf('  %s\n', repmat('-', 1, 64));
        fprintf('  %-28s  %10.3f  %10.3f  %+9.1f%%\n', ...
            sprintf('Peak (%s)', mod.ylabel), m2d.peak_val, m3d.peak_val, 100*(peak_ratio-1));
        fprintf('  %-28s  %10.2f  %10.2f  %+9.2f mm\n', ...
            'Focal distance (mm)', m2d.focal_dist_mm, m3d.focal_dist_mm, focal_diff_mm);
        fprintf('  %-28s  %10.2f  %10.2f  %+9.2f mm\n', ...
            'FWHM axial (mm)',     m2d.fwhm_mm,       m3d.fwhm_mm,       fwhm_diff_mm);
        fprintf('  %-28s  %10.4f\n', 'Pearson r (2D vs 3D)', r);

        % Accumulate table columns
        pfx       = mod.name(1:3);
        row_data  = [row_data, {m2d.peak_val, m3d.peak_val, peak_ratio, ...
                                m2d.focal_dist_mm, m3d.focal_dist_mm, focal_diff_mm, ...
                                m2d.fwhm_mm, m3d.fwhm_mm, fwhm_diff_mm, r}]; %#ok<AGROW>
        row_vars  = [row_vars, { ...
            [pfx '_peak_2D'],  [pfx '_peak_3D'],  [pfx '_peak_ratio'], ...
            [pfx '_focal_2D'], [pfx '_focal_3D'], [pfx '_focal_diff'], ...
            [pfx '_fwhm_2D'],  [pfx '_fwhm_3D'],  [pfx '_fwhm_diff'], ...
            [pfx '_pearson_r']}]; %#ok<AGROW>

        % ------------------------------------------------------------------
        % Figures: absolute and peak-normalised profiles
        % ------------------------------------------------------------------
        fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 22 9], ...
                     'Name', sprintf('%s — %s', bm.label, mod.name));

        subplot(1, 2, 1);
        plot(z_2d, prof2d, 'b-',  'LineWidth', 1.5, 'DisplayName', '2D axisymmetric'); hold on;
        plot(z_3d, prof3d, 'r--', 'LineWidth', 1.5, 'DisplayName', '3D full');
        xlabel('Distance from transducer face (mm)');
        ylabel(mod.ylabel, 'Interpreter', 'tex');
        title(sprintf('%s — %s (absolute)', bm.label, mod.name), 'Interpreter', 'none');
        legend('Location', 'best');
        xlim([0, max(z_2d(end), z_3d(end))]);
        grid on; box on;

        subplot(1, 2, 2);
        plot(z_2d, prof2d / m2d.peak_val, 'b-',  'LineWidth', 1.5, 'DisplayName', '2D axisymmetric'); hold on;
        plot(z_3d, prof3d / m3d.peak_val, 'r--', 'LineWidth', 1.5, 'DisplayName', '3D full');
        xline(m2d.focal_dist_mm, 'b:', 'LineWidth', 1);
        xline(m3d.focal_dist_mm, 'r:', 'LineWidth', 1);
        xlabel('Distance from transducer face (mm)');
        ylabel('Normalised');
        title(sprintf('Peak-normalised  (r = %.3f)', r), 'Interpreter', 'none');
        legend('Location', 'best');
        xlim([0, max(z_2d(end), z_3d(end))]);
        ylim([0 1.1]);
        grid on; box on;

        figname = fullfile(pn.figures, sprintf('compare_2D_3D_%s_%s', sub_str, mod.name));
        saveas(fig, [figname, '.png']);
        saveas(fig, [figname, '.eps'], 'epsc');
        fprintf('\n  Figure saved: %s.png/.eps\n', figname);
        close(fig);
    end

    if numel(row_data) > 2
        results_table = [results_table; cell2table(row_data, 'VariableNames', row_vars)]; %#ok<AGROW>
    end
end

% ------------------------------------------------------------------
% Save summary table
% ------------------------------------------------------------------
if ~isempty(results_table)
    csv_out = fullfile(pn.figures, 'compare_2D_3D_metrics.csv');
    writetable(results_table, csv_out);
    fprintf('\nSummary saved: %s\n', csv_out);
end

fprintf('\nDone.\n');

% -------------------------------------------------------------------------
% Local function
% -------------------------------------------------------------------------
function m = profile_metrics(prof, z_mm)
    [pk, pk_idx] = max(prof);
    m.peak_val      = pk;
    m.focal_dist_mm = z_mm(pk_idx);
    half_pk = pk / 2;
    above   = prof >= half_pk;
    tr      = diff(above);
    ri      = find(tr ==  1, 1, 'first');
    fi      = find(tr == -1, 1, 'last');
    if ~isempty(ri) && ~isempty(fi)
        m.fwhm_mm = z_mm(fi) - z_mm(ri);
    else
        m.fwhm_mm = NaN;
    end
end
