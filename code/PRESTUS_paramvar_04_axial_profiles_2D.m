% PRESTUS_PARAMVAR_04_AXIAL_PROFILES_2D
%
% Extracts and compares focal-axis intensity AND thermal profiles for all 2D
% axisymmetric parameter-variation runs in sub-001 (water/brain) and
% sub-002 (skull/brain), at the primary 0.5 mm grid resolution.
%
% For each tissue × property combination, the script plots the axial profile
% for the bottom (_1) and top (_2) sweep levels superimposed, showing how
% each tissue parameter affects the focal depth and peak values.
%
% Acoustic modality  (always present):
%   NII suffix: _Ipa.nii.gz
%   Metrics: peak ISPPA (W/cm²), focal distance (mm), FWHM (mm)
%
% Thermal modality   (present once thermal sims complete):
%   NII suffix: _heatrise.nii.gz   — peak temperature rise (°C above baseline)
%   Metrics: peak ΔT (°C), location of peak ΔT (mm), FWHM of thermal profile (mm)
%   Note: If thermal NIfTIs are absent the thermal panels are skipped with a
%         warning.  Re-run simulations with run_heating_sims:1 and ensure
%         nifti_thermal completes (check log for "400/400 pulse trains").
%
% Output figures are saved to figures/paramvar/axial_profiles/.
%
% See also: PRESTUS_PARAMVAR_01_BENCHMARKS_SIM

clear; close all; clc;

% -------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
addpath('/project/2425122.01/v05/tools/nifti_toolbox');
if isfolder('/Volumes/2425122.01/v05/tools/nifti_toolbox')
    addpath('/Volumes/2425122.01/v05/tools/nifti_toolbox');
end

pn.figures = fullfile(rootpath, 'figures', 'paramvar', 'axial_profiles');
if ~isfolder(pn.figures); mkdir(pn.figures); end

dir_2d = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a');

% -------------------------------------------------------------------------
% Configuration
% -------------------------------------------------------------------------
voxel_mm  = 0.5;
trans_vox = 1;    % transducer face at z-voxel 1 (config: trans_pos z = 1)

subjects = struct( ...
    'id',              {1,             2}, ...
    'label',           {'water/brain', 'skull/brain'}, ...
    'tissue_labels',   {{'brain'},     {'brain', 'skull'}}, ...
    'tissue_suffixes', {{'_b'},        {'_b', '_s'}});

property_labels   = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                     'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};
property_symbols  = {'c', '\rho', '\alpha_0', '\kappa', 'c_p', '\omega', '\beta'};

level_colors = {[0.2 0.5 0.9], [0.9 0.3 0.2]};   % blue = bottom, red = top
level_labels = {'bottom (min)', 'top (max)'};

% Modalities: each specifies the NII suffix and axis label
modalities = struct( ...
    'name',    {'acoustic',           'thermal'}, ...
    'suffix',  {'_Ipa.nii.gz',        '_heatrise.nii.gz'}, ...
    'ylabel',  {'Intensity (W/cm²)',  '\DeltaT (°C)'}, ...
    'metric',  {'peak\_Isppa',        'peak\_dT'});

% -------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------
for s = 1:numel(subjects)
    subj    = subjects(s);
    sub_id  = subj.id;
    sub_str = sprintf('sub-%03d', sub_id);
    nii_dir = fullfile(dir_2d, sub_str, 'nii');

    fprintf('\n============================================================\n');
    fprintf('  %s  (%s)\n', sub_str, subj.label);
    fprintf('============================================================\n');

    % Determine grid parameters from first available intensity NII
    ref_files = dir(fullfile(nii_dir, sprintf('%s_*_Ipa.nii.gz', sub_str)));
    if isempty(ref_files)
        warning('No intensity NII files found for %s — skipping.', sub_str);
        continue;
    end
    ref_img = load_nii(fullfile(nii_dir, ref_files(1).name));
    ref_vol = double(ref_img.img);
    % Phantoms are true 2-D arrays [x (lateral), z (axial)] — no singleton
    % y-dimension (see createPhantom.m array-layout header) — so the axial
    % axis is dim 2, not dim 3.
    cx   = round(size(ref_vol, 1) / 2);
    nz   = size(ref_vol, 2);
    z_mm = ((1:nz) - trans_vox) * voxel_mm;

    % Check whether thermal NIfTIs exist for this subject
    thermal_files = dir(fullfile(nii_dir, sprintf('%s_*_heatrise.nii.gz', sub_str)));
    has_thermal   = ~isempty(thermal_files);
    if ~has_thermal
        warning(['No heatrise NII files found for %s.\n' ...
                 'Thermal panels will be skipped.\n' ...
                 'Thermal sims may still be running — check log for 400/400 pulse trains,\n' ...
                 'then re-run this script once nifti_thermal has completed.'], sub_str);
    end

    active_modalities = modalities(1);   % acoustic always present
    if has_thermal
        active_modalities = modalities;  % add thermal
    end
    n_mod = numel(active_modalities);

    % Collect summary rows per modality
    summaries = cell(1, n_mod);
    for mo = 1:n_mod
        summaries{mo} = struct('label', {}, 'level', {}, ...
                               'peak_val', {}, 'focal_dist_mm', {}, 'fwhm_mm', {});
    end

    % ------------------------------------------------------------------
    % Per-tissue profile panels
    % ------------------------------------------------------------------
    for t = 1:numel(subj.tissue_labels)
        tissue        = subj.tissue_labels{t};
        tissue_suffix = subj.tissue_suffixes{t};
        n_props       = numel(property_labels);

        for mo = 1:n_mod
            mod       = active_modalities(mo);
            fig = figure('Color', 'w', 'Units', 'centimeters', ...
                         'Position', [1 1 28 20], ...
                         'Name', sprintf('%s — %s — %s tissue', sub_str, mod.name, tissue));

            for p = 1:n_props
                prop_suffix = property_suffixes{p};
                prop_sym    = property_symbols{p};

                ax = subplot(ceil(n_props/4), 4, p);
                hold(ax, 'on');

                for lv = 1:2
                    nii_path = fullfile(nii_dir, sprintf('%s_phantom_T1w%s%s%d%s', ...
                        sub_str, tissue_suffix, prop_suffix, lv, mod.suffix));

                    if ~isfile(nii_path)
                        warning('File not found: %s', nii_path);
                        continue;
                    end

                    nii  = load_nii(nii_path);
                    vol  = double(nii.img);
                    prof = vol(cx, :)';
                    prof(prof < 0) = 0;

                    plot(ax, z_mm, prof, 'Color', level_colors{lv}, ...
                         'LineWidth', 1.5, 'DisplayName', level_labels{lv});

                    m = profile_metrics(prof, z_mm);

                    row.label        = sprintf('%s%s%s%d', sub_str, tissue_suffix, prop_suffix, lv);
                    row.level        = lv;
                    row.peak_val     = m.peak_val;
                    row.focal_dist_mm = m.focal_dist_mm;
                    row.fwhm_mm      = m.fwhm_mm;
                    summaries{mo}(end+1) = row;

                    fprintf('  [%s] %-35s  peak=%7.3f  focal=%5.1f mm  FWHM=%5.1f mm\n', ...
                        mod.name, row.label, m.peak_val, m.focal_dist_mm, m.fwhm_mm);
                end

                xlabel(ax, 'Distance from transducer (mm)', 'FontSize', 8);
                ylabel(ax, mod.ylabel, 'Interpreter', 'tex', 'FontSize', 8);
                title(ax, prop_sym, 'Interpreter', 'tex', 'FontSize', 10);
                grid(ax, 'on'); box(ax, 'on');
                xlim(ax, [0, z_mm(end)]);
                if p == 1
                    legend(ax, 'Location', 'best', 'FontSize', 7);
                end
            end

            sgtitle(sprintf('%s — %s — %s tissue — axial profiles', ...
                            sub_str, mod.name, tissue), 'Interpreter', 'none');

            figname = fullfile(pn.figures, sprintf('%s_%s_%s_axial_profiles', ...
                sub_str, mod.name, tissue));
            saveas(fig, [figname, '.png']);
            saveas(fig, [figname, '.eps'], 'epsc');
            fprintf('\n  Saved: %s.png/.eps\n', figname);
            close(fig);
        end
    end

    % ------------------------------------------------------------------
    % Combined summary figure (acoustic + thermal side by side if available)
    % ------------------------------------------------------------------
    for mo = 1:n_mod
        sm  = summaries{mo};
        mod = active_modalities(mo);
        if isempty(sm); continue; end

        n          = numel(sm);
        labels_all = {sm.label};
        pk_all     = [sm.peak_val];
        fd_all     = [sm.focal_dist_mm];
        lv_all     = [sm.level];
        bar_colors = level_colors(lv_all);

        fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [1 1 32 10], ...
                      'Name', sprintf('%s — %s summary', sub_str, mod.name));

        subplot(1, 3, 1);
        bh = bar(pk_all);
        bh.FaceColor = 'flat';
        for k = 1:n; bh.CData(k,:) = bar_colors{k}; end
        set(gca, 'XTick', 1:n, ...
            'XTickLabel', strrep(labels_all, sub_str, ''), ...
            'XTickLabelRotation', 45, 'FontSize', 7);
        ylabel(mod.ylabel, 'Interpreter', 'tex');
        title(sprintf('Peak (%s)', mod.name));
        grid on;

        subplot(1, 3, 2);
        bh2 = bar(fd_all);
        bh2.FaceColor = 'flat';
        for k = 1:n; bh2.CData(k,:) = bar_colors{k}; end
        set(gca, 'XTick', 1:n, ...
            'XTickLabel', strrep(labels_all, sub_str, ''), ...
            'XTickLabelRotation', 45, 'FontSize', 7);
        ylabel('Focal distance (mm)');
        title('Depth of peak');
        grid on;

        subplot(1, 3, 3);
        fwhm_all = [sm.fwhm_mm];
        bh3 = bar(fwhm_all);
        bh3.FaceColor = 'flat';
        for k = 1:n; bh3.CData(k,:) = bar_colors{k}; end
        set(gca, 'XTick', 1:n, ...
            'XTickLabel', strrep(labels_all, sub_str, ''), ...
            'XTickLabelRotation', 45, 'FontSize', 7);
        ylabel('FWHM (mm)');
        title('Axial FWHM');
        grid on;

        sgtitle(sprintf('%s — %s — parameter sensitivity summary', sub_str, mod.name), ...
                'Interpreter', 'none');
        figname2 = fullfile(pn.figures, sprintf('%s_%s_summary_metrics', sub_str, mod.name));
        saveas(fig2, [figname2, '.png']);
        saveas(fig2, [figname2, '.eps'], 'epsc');
        fprintf('  Saved: %s.png/.eps\n', figname2);
        close(fig2);

        % Save CSV
        T       = struct2table(sm);
        csv_out = fullfile(pn.figures, sprintf('%s_%s_axial_profile_metrics.csv', sub_str, mod.name));
        writetable(T, csv_out);
        fprintf('  Saved: %s\n', csv_out);
    end
end

fprintf('\nDone.\n');

% -------------------------------------------------------------------------
% Local function
% -------------------------------------------------------------------------
function m = profile_metrics(prof, z_mm)
    [pk, pk_idx] = max(prof);
    m.peak_val      = pk;
    m.focal_dist_mm = z_mm(pk_idx);
    half_pk  = pk / 2;
    above    = prof >= half_pk;
    tr       = diff(above);
    ri       = find(tr ==  1, 1, 'first');
    fi       = find(tr == -1, 1, 'last');
    if ~isempty(ri) && ~isempty(fi)
        m.fwhm_mm = z_mm(fi) - z_mm(ri);
    else
        m.fwhm_mm = NaN;
    end
end
