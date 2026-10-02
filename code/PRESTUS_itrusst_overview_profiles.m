% PRESTUS_ITRUSST_OVERVIEW_PROFILES
%
% Loads all available itrusst simulation NIfTIs and produces two figures
% per analysis folder (one folder shown at a time):
%
%   Band plot
%     Axial intensity and heatrise profiles. The band spans min to max
%     across all paramvar variants; the centre line is the mean.
%     Both subjects are overlaid (solid = sub-001, dashed = sub-002).
%     Tissue transitions are marked with vertical dashed lines.
%
%   Butterfly plot
%     Same layout but every individual variant profile is drawn as a thin
%     semi-transparent line; the mean is overlaid as a thick line.
%
% Subjects covered: sub-001 (water/brain), sub-002 (skull/brain).
% Extend SUB_DEFS below to include other subjects.
%
% Output: figures/paramvar/overview_profiles/{folder}_{band|butterfly}.png

%% ── Paths ───────────────────────────────────────────────────────────────────

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'));
rootpath = pwd;

pn.tussim   = fullfile(rootpath, 'data', 'tussim');
pn.simnibs  = fullfile(rootpath, 'data', 'simnibs');
pn.figures  = fullfile(rootpath, 'figures', 'paramvar', 'overview_profiles');
if ~isfolder(pn.figures); mkdir(pn.figures); end

%% ── Subject definitions ─────────────────────────────────────────────────────

SUB_DEFS = struct( ...
    'id',        {1,             2}, ...
    'label',     {'water/brain', 'skull/brain'}, ...
    'linestyle', {'-',           '--'});

%% ── Output folder definitions ───────────────────────────────────────────────
% Each row: {folder_name, display_label, RGB_color}

FOLDERS = { ...
    'itrusst_protocol1a',      '0.5 mm param sweep',   [0.18 0.45 0.78]; ...
};

%% ── Modalities ──────────────────────────────────────────────────────────────

MODS = struct( ...
    'suffix', {'_Ipa.nii.gz',       '_heatrise.nii.gz'}, ...
    'ylabel', {'I_{SPPA} (W/cm²)', '\DeltaT (°C)'}, ...
    'name',   {'intensity',         'heatrise'});

%% ── Tissue label colours for transition annotations ─────────────────────────

TISSUE_COLORS = containers.Map( ...
    {0, 1, 4, 7, 8}, ...
    {[0.75 0.90 0.95], [0.70 0.85 0.70], [0.95 0.82 0.60], ...
     [0.95 0.82 0.60], [0.92 0.90 0.70]});
TISSUE_NAMES  = containers.Map( ...
    {0, 1, 4, 7, 8}, ...
    {'water', 'brain', 'skull', 'cortical', 'trabecular'});

%% ── Collect data for all subjects ───────────────────────────────────────────

n_subs    = numel(SUB_DEFS);
n_folders = size(FOLDERS, 1);
n_mods    = numel(MODS);

% all_data{s,f,m}  — struct with .profiles [n_var×nz] and .z_mm [1×nz]
all_data  = cell(n_subs, n_folders, n_mods);
all_trans = cell(n_subs, 1);

for s = 1:n_subs
    sub_id  = SUB_DEFS(s).id;
    sub_str = sprintf('sub-%03d', sub_id);
    fprintf('\n[%s]  %s\n', sub_str, SUB_DEFS(s).label);

    phantom_path = fullfile(pn.simnibs, sprintf('m2m_%s', sub_str), 'final_tissues.nii.gz');
    all_trans{s} = load_tissue_transitions(phantom_path);

    for f = 1:n_folders
        folder_name = FOLDERS{f, 1};
        nii_dir     = fullfile(pn.tussim, folder_name, sub_str, 'nii');
        if ~isfolder(nii_dir)
            fprintf('  [skip] %s — folder not found\n', folder_name);
            continue;
        end

        for m = 1:n_mods
            pattern = sprintf('%s_*%s', sub_str, MODS(m).suffix);
            files   = dir(fullfile(nii_dir, pattern));

            base_pat = sprintf('%s_phantom_T1w%s', sub_str, MODS(m).suffix);
            base_f   = fullfile(nii_dir, base_pat);
            if isfile(base_f) && isempty(files)
                files = dir(base_f);
            end

            if isempty(files); continue; end

            profiles = [];
            z_mm     = [];
            for k = 1:numel(files)
                fpath = fullfile(nii_dir, files(k).name);
                [prof, zmm] = extract_axial_profile(fpath);
                if isempty(prof); continue; end
                if isempty(profiles)
                    profiles = zeros(numel(files), numel(prof));
                    z_mm     = zmm;
                end
                profiles(k, :) = prof(:)';
            end

            valid    = any(profiles ~= 0, 2);
            profiles = profiles(valid, :);

            if ~isempty(profiles)
                all_data{s, f, m}.profiles = profiles;
                all_data{s, f, m}.z_mm     = z_mm;
                fprintf('  %-30s  %-10s  %d variants\n', ...
                    folder_name, MODS(m).name, size(profiles, 1));
            end
        end
    end
end

%% ── Per-folder figures (both subjects overlaid) ──────────────────────────────
%
% One band figure + one butterfly figure per analysis folder per modality.
% Color   = subject  (sub-001: black, sub-002: red)
% Linestyle = solid for mean; shaded band for min–max range.

SUB_COLORS = [0.00 0.00 0.00; ...   % sub-001  black
              0.85 0.10 0.10];       % sub-002  red

FONT_SIZE = 13;

for f = 1:n_folders
    folder_name  = FOLDERS{f, 1};
    folder_label = FOLDERS{f, 2};

    % Skip if no subject has data in this folder for any modality
    has_data = false;
    for s = 1:n_subs
        for m = 1:n_mods
            if ~isempty(all_data{s, f, m}); has_data = true; break; end
        end
        if has_data; break; end
    end
    if ~has_data
        fprintf('  [skip figure] %s — no data for any subject\n', folder_name);
        continue;
    end

    for m = 1:n_mods
        mod_has_data = false;
        for s = 1:n_subs
            if ~isempty(all_data{s, f, m}); mod_has_data = true; break; end
        end
        if ~mod_has_data; continue; end

        fig_band = figure('Color', 'w', 'Units', 'centimeters', ...
            'Position', [2 2 24 10], ...
            'Name', sprintf('%s — %s band', folder_name, MODS(m).name));
        set(fig_band, 'Renderer', 'painters');
        set(0, 'DefaultFigureRenderer', 'painters');

        fig_butt = figure('Color', 'w', 'Units', 'centimeters', ...
            'Position', [2 2 24 10], ...
            'Name', sprintf('%s — %s butterfly', folder_name, MODS(m).name));
        set(fig_butt, 'Renderer', 'painters');

        ax_band = axes('Parent', fig_band);
        ax_butt = axes('Parent', fig_butt);
        hold(ax_band, 'on'); hold(ax_butt, 'on');

        x_max = 0;

        % Pass 1: fills and butterfly individual lines (below mean lines)
        for s = 1:n_subs
            if isempty(all_data{s, f, m}); continue; end
            d    = all_data{s, f, m};
            clr  = SUB_COLORS(s, :);
            lo   = min(d.profiles, [], 1);
            hi   = max(d.profiles, [], 1);
            n_var = size(d.profiles, 1);

            % Pre-blend colors with white (painters renderer has no alpha support)
            band_clr = 1 - (1 - clr) * 0.20;
            alpha_val = max(0.12, min(0.35, 1 / sqrt(n_var)));
            butt_clr  = 1 - (1 - clr) * alpha_val;

            fill(ax_band, [d.z_mm, fliplr(d.z_mm)], [lo, fliplr(hi)], band_clr, ...
                'EdgeColor', 'none');
            for k = 1:n_var
                plot(ax_butt, d.z_mm, d.profiles(k, :), '-', ...
                    'Color', butt_clr, 'LineWidth', 0.8);
            end

            x_max = max(x_max, d.z_mm(end));
        end

        % Pass 2: mean lines on top of all fills
        for s = 1:n_subs
            if isempty(all_data{s, f, m}); continue; end
            d   = all_data{s, f, m};
            clr = SUB_COLORS(s, :);
            mn  = mean(d.profiles, 1);
            plot(ax_band, d.z_mm, mn, '-', 'Color', clr, 'LineWidth', 1.8);
            plot(ax_butt, d.z_mm, mn, '-', 'Color', clr, 'LineWidth', 2.2);
        end

        if x_max == 0
            close(fig_band); close(fig_butt);
            continue;
        end

        for ax = [ax_band, ax_butt]
            xlim(ax, [0, x_max]);
            for s = 1:n_subs
                add_tissue_transitions(ax, all_trans{s}, TISSUE_COLORS, TISSUE_NAMES);
            end
            xlabel(ax, '');
            yl = ylabel(ax, MODS(m).ylabel, 'Interpreter', 'tex', 'FontSize', FONT_SIZE);
            set(yl, 'Rotation', 90);
            ax.FontSize = FONT_SIZE;
            grid(ax, 'off'); box(ax, 'on');
        end

        sgtitle(fig_band, sprintf('%s — mean ± range', folder_label), ...
            'Interpreter', 'none', 'FontSize', FONT_SIZE);
        sgtitle(fig_butt, sprintf('%s — butterfly: all paramvar profiles', folder_label), ...
            'Interpreter', 'none', 'FontSize', FONT_SIZE);

        save_figure(fig_band, fullfile(pn.figures, sprintf('%s_band_%s',      folder_name, MODS(m).name)));
        save_figure(fig_butt, fullfile(pn.figures, sprintf('%s_butterfly_%s', folder_name, MODS(m).name)));
        fprintf('  Saved figures for %s — %s\n', folder_name, MODS(m).name);
        close(fig_band); close(fig_butt);
    end
end

fprintf('\nDone.\n');

%% ═══════════════════════════════════════════════════════════════════════════
%% Local functions
%% ═══════════════════════════════════════════════════════════════════════════

function save_figure(fig, base_path)
% SAVE_FIGURE  Export figure as PNG (300 dpi) and EPS (vector, Illustrator-ready).
% base_path should have no extension — .png and .eps are appended.
%
% Painters renderer is set on the figure before calling this function.
% EPS color (epsc) is fully editable in Illustrator; no transparency issues.

    % PNG — raster for quick viewing
    print(fig, [base_path '.png'], '-dpng', '-r300');

    % EPS — vector, editable in Illustrator
    saveas(fig, [base_path '.eps'], 'epsc');
end

% ─────────────────────────────────────────────────────────────────────────────

function [prof, z_mm] = extract_axial_profile(nii_path)
% EXTRACT_AXIAL_PROFILE  Read a NIfTI and return the beam-axis intensity
% vector and the corresponding depth axis in mm.
%
% Handles 2D [nx × nz] and 3D [nx × ny × nz] volumes. The beam axis is at
% the lateral centre (dim1 = dim2 = round(n/2)). The z_mm axis starts at 0
% at the transducer face; trans_vox is inferred from voxel size.

    prof  = [];
    z_mm  = [];

    try
        info = niftiinfo(nii_path);
        vol  = double(niftiread(nii_path));
    catch err
        warning('Failed to read %s: %s', nii_path, err.message);
        return;
    end

    vol(vol < 0) = 0;
    sz  = size(vol);
    res = info.PixelDimensions(end);   % axial voxel size [mm] (last dim)

    % trans_vox: first grid voxel occupied by the transducer face
    % 0.5 mm grid → trans_pos z = 1; 0.25 mm grid → trans_pos z = 4
    if res >= 0.4
        trans_vox = 1;
    else
        trans_vox = 4;
    end

    switch numel(sz)
        case 2
            % [nx × nz]: lateral × axial
            cx   = round(sz(1) / 2);
            prof = vol(cx, :);
            nz   = sz(2);
        case 3
            % [nx × ny × nz]: x × y × axial (radially expanded)
            cx   = round(sz(1) / 2);
            cy   = round(sz(2) / 2);
            prof = squeeze(vol(cx, cy, :));
            nz   = sz(3);
        otherwise
            warning('Unexpected NIfTI ndim=%d in %s', numel(sz), nii_path);
            return;
    end

    z_mm = ((1:nz) - trans_vox) * res;
    prof = prof(:)';
end

% ─────────────────────────────────────────────────────────────────────────────

function trans = load_tissue_transitions(phantom_path)
% LOAD_TISSUE_TRANSITIONS  Read a 2D phantom NIfTI and return a struct array
% of tissue transitions along the beam axis.
%
% Each element: trans(k).z_mm, trans(k).from_label, trans(k).to_label

    trans = struct('z_mm', {}, 'from_label', {}, 'to_label', {});

    if ~isfile(phantom_path)
        warning('Phantom not found: %s', phantom_path);
        return;
    end

    try
        info = niftiinfo(phantom_path);
        vol  = round(double(niftiread(phantom_path)));
    catch
        return;
    end

    vol = squeeze(vol);
    sz  = size(vol);
    res = info.PixelDimensions(1);

    if numel(sz) == 2
        cx     = round(sz(1) / 2);
        labels = vol(cx, :);
    else
        return;
    end

    trans_vox = 1;
    changes   = find(diff(labels) ~= 0);
    for k = 1:numel(changes)
        t.z_mm       = (changes(k) - trans_vox + 0.5) * res;
        t.from_label = labels(changes(k));
        t.to_label   = labels(changes(k) + 1);
        trans(end+1) = t; %#ok<AGROW>
    end
end

% ─────────────────────────────────────────────────────────────────────────────

function add_tissue_transitions(ax, trans, tissue_colors, tissue_names)
% ADD_TISSUE_TRANSITIONS  Draw vertical dashed lines and tissue-region labels.
% Called after xlim is set; uses normalised axes coordinates for label
% placement so it is robust to y-scale changes.

    if isempty(trans); return; end

    xlims = xlim(ax);

    for k = 1:numel(trans)
        z = trans(k).z_mm;
        if z < xlims(1) || z > xlims(2); continue; end

        xline(ax, z, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, ...
            'HandleVisibility', 'off');

        if tissue_names.isKey(trans(k).to_label)
            text(ax, z + 0.5 * diff(xlims) / 100, ...
                ax.YLim(1) + 0.90 * diff(ax.YLim), ...
                tissue_names(trans(k).to_label), ...
                'FontSize', 7, 'Color', [0.35 0.35 0.35], ...
                'Rotation', 90, 'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'left', 'HandleVisibility', 'off');
        end
    end
end
