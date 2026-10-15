% PAPER_FIG23_AXIAL_PROFILES
%
% Paper Figs. 2 and 3 (second panel): on-axis I_SPPA (Fig. 2) and heat
% rise (Fig. 3) profiles, drawn vertically (depth downwards, transducer at
% the top) to sit next to the benchmark schematic. Mean (line) and
% min-max range (band) across all one-at-a-time parameter variations of
% itrusst_protocol1a, for the benchmark without skull (sub-001, black)
% and with skull (sub-002, red). A strip on the right shows the sub-002
% tissue layers along the axis (water / skull / brain).
%
% Profile extraction as in PRESTUS_itrusst_overview_profiles.m.
%
% Output: figures/paper/fig2/fig2_axial_profile_intensity.{pdf,png}
%         figures/paper/fig3/fig3_axial_profile_heatrise.{pdf,png}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn.tussim  = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a');
pn.simnibs = fullfile(rootpath, 'data', 'simnibs');
pn.figures = fullfile(rootpath, 'figures', 'paper');

%% settings

depth_lim = [0 100]; % mm from the transducer
mods = struct( ...
    'suffix', {'_Ipa.nii.gz',            '_heatrise.nii.gz'}, ...
    'xlabel', {'I_{SPPA} (W/cm^2)',      'Temperature rise (°C)'}, ...
    'fig',    {'fig2',                   'fig3'}, ...
    'name',   {'fig2_axial_profile_intensity', 'fig3_axial_profile_heatrise'});
subjects = [1 2];

S = paper_style();
sub_colors = {S.noskull, S.withskull};

%% tissue layers along the axis (sub-002)

tissue_info = niftiinfo(fullfile(pn.simnibs, 'm2m_sub-002', 'final_tissues.nii.gz'));
tissue = round(double(niftiread(tissue_info)));
tissue_axis = tissue(round(size(tissue, 1) / 2), :);
tissue_z = ((1:numel(tissue_axis)) - 1) * tissue_info.PixelDimensions(end);
transitions = tissue_z(find(diff(tissue_axis) ~= 0)) + tissue_info.PixelDimensions(end) / 2;

%% plot

for m = 1:numel(mods)
    h = paper_figure(S.profile_fig);
    ax = axes(h, 'Units', 'centimeters', 'Position', S.profile_axes);
    hold(ax, 'on');

    x_max = 0;
    profiles = cell(1, numel(subjects));
    for s = 1:numel(subjects)
        sub_str = sprintf('sub-%03d', subjects(s));
        files = dir(fullfile(pn.tussim, sub_str, 'nii', [sub_str '_phantom_T1w*' mods(m).suffix]));
        P = [];
        for k = 1:numel(files)
            [prof, z_mm] = extract_axial_profile(fullfile(files(k).folder, files(k).name));
            P(end+1, :) = prof; %#ok<AGROW>
        end
        if isempty(P)
            warning('%s: no %s profiles found.', sub_str, mods(m).suffix);
            continue;
        end
        fprintf('%s %s: %d profiles\n', sub_str, mods(m).suffix, size(P, 1));
        profiles{s} = struct('P', P, 'z', z_mm);
        band_clr = 1 - (1 - sub_colors{s}) * 0.2; % pre-blended (no alpha in painters)
        fill(ax, [min(P, [], 1), fliplr(max(P, [], 1))], [z_mm, fliplr(z_mm)], band_clr, 'EdgeColor', 'none');
        x_max = max(x_max, max(P(:, z_mm <= depth_lim(2)), [], 'all'));
    end
    for s = 1:numel(subjects)
        if isempty(profiles{s}); continue; end
        plot(ax, mean(profiles{s}.P, 1), profiles{s}.z, '-', 'Color', sub_colors{s}, 'LineWidth', S.linewidth);
    end

    x_lim = [0, ceil(x_max / 10^floor(log10(x_max)) ) * 10^floor(log10(x_max))];
    for z = transitions(transitions > depth_lim(1) & transitions < depth_lim(2))
        plot(ax, x_lim, [z z], '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.5);
    end
    set(ax, 'YDir', 'reverse', 'YLim', depth_lim, 'XLim', x_lim, 'YTick', [], ...
        'YColor', 'none', 'Layer', 'top');
    title(ax, mods(m).xlabel, 'Interpreter', 'tex');

    % tissue strip
    axs = axes(h, 'Units', 'centimeters', 'Position', ...
        [S.profile_axes(1) + S.profile_axes(3) + 0.1, S.profile_axes(2), 0.18, S.profile_axes(4)]);
    tissue_clr = containers.Map({0, 1, 4}, {S.water, S.brain, S.skull});
    seg = [0, find(diff(tissue_axis) ~= 0), numel(tissue_axis)];
    hold(axs, 'on');
    for k = 1:numel(seg) - 1
        z0 = tissue_z(seg(k) + 1); z1 = tissue_z(seg(k + 1)) + tissue_info.PixelDimensions(end);
        lab = tissue_axis(seg(k) + 1);
        if ~tissue_clr.isKey(lab); continue; end
        patch(axs, [0 1 1 0], [z0 z0 z1 z1], tissue_clr(lab), 'EdgeColor', 'none');
    end
    set(axs, 'YDir', 'reverse', 'YLim', depth_lim, 'XLim', [0 1], 'Visible', 'off');

    paper_save(h, fullfile(pn.figures, mods(m).fig), mods(m).name);
    close(h);
end

%% local functions

function [prof, z_mm] = extract_axial_profile(nii_path)
% beam-axis profile of a 2D [lateral x axial] NIfTI; z = 0 at the
% transducer (first grid voxel), as in PRESTUS_itrusst_overview_profiles.m
info = niftiinfo(nii_path);
vol  = double(niftiread(info));
vol(vol < 0) = 0;
res  = info.PixelDimensions(end);
prof = vol(round(size(vol, 1) / 2), :);
z_mm = ((1:numel(prof)) - 1) * res;
end
