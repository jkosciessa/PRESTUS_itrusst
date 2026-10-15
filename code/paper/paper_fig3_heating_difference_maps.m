% PAPER_FIG3_HEATING_DIFFERENCE_MAPS
%
% Paper Fig. 3 (right): maps of the heat-rise difference between the top
% and bottom level of a skull property (sub-002, itrusst_protocol1a), with
% the skull outline. Shown for skull absorption fraction (Fig. 3a) and
% skull attenuation (Fig. 3b). The colour range is symmetric around zero
% (black: cooler, red: warmer with the top level).
%
% Transducer, beam cone and frames are added in Illustrator.
%
% Output: figures/paper/fig3/fig3_diffmap_<property>.{pdf,png}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn.data    = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a', 'sub-002');
pn.simnibs = fullfile(rootpath, 'data', 'simnibs');
pn.figures = fullfile(rootpath, 'figures', 'paper', 'fig3');

%% settings

maps = struct( ...
    'affix', {'_s_abs_',                    '_s_a0_'}, ...
    'title', {'skull absorption fraction',  'skull attenuation'}, ...
    'name',  {'fig3_diffmap_skull_absorption', 'fig3_diffmap_skull_attenuation'});
depth_lim   = [0 90];   % mm from the transducer
lateral_lim = [-35 35]; % mm from the axis

S = paper_style();
cmap = [linspace(0, 1, 128)', linspace(0, 1, 128)', linspace(0, 1, 128)'; ...    % black -> white
        ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)'];                % white -> red

%% skull outline

tissue_info = niftiinfo(fullfile(pn.simnibs, 'm2m_sub-002', 'final_tissues.nii.gz'));
skull = double(round(niftiread(tissue_info)) == 4);

%% plot

for m = 1:numel(maps)
    lvl = cell(1, 2);
    for l = 1:2
        f = fullfile(pn.data, 'nii', sprintf('sub-002_phantom_T1w%s%d_heatrise.nii.gz', maps(m).affix, l));
        lvl{l} = double(niftiread(f));
    end
    d = lvl{2} - lvl{1}; % top - bottom
    info = niftiinfo(f);
    res = info.PixelDimensions(1);
    x_mm = ((1:size(d, 1)) - round(size(d, 1) / 2)) * res;
    z_mm = ((1:size(d, 2)) - 1) * res;
    c_lim = max(abs(d(:)));
    c_lim = ceil(c_lim * 10) / 10;
    fprintf('%s: top - bottom in [%.3f, %.3f] °C\n', maps(m).title, min(d(:)), max(d(:)));

    h = paper_figure([5.2 6.6]);
    ax = axes(h, 'Units', 'centimeters', 'Position', [0.3 0.5 3.6 3.6 * diff(depth_lim) / diff(lateral_lim)]);
    imagesc(ax, x_mm, z_mm, d', [-c_lim c_lim]);
    hold(ax, 'on');
    contour(ax, x_mm, z_mm, skull', [0.5 0.5], 'Color', [0.2 0.2 0.2], 'LineWidth', 0.5);
    colormap(ax, cmap);
    axis(ax, 'image');
    set(ax, 'XLim', lateral_lim, 'YLim', depth_lim, 'YDir', 'reverse', 'Visible', 'off');
    title(ax, maps(m).title, 'Visible', 'on');

    cb = colorbar(ax, 'Units', 'centimeters', 'Position', [4.2 2.0 0.25 2.5]);
    cb.Ticks = [-c_lim 0 c_lim];
    cb.Label.String = {'Heating difference (°C)', '(top - bottom)'};
    cb.TickDirection = 'out';
    cb.LineWidth = 0.5;
    set(ax, 'Position', [0.3 0.5 3.6 3.6 * diff(depth_lim) / diff(lateral_lim)]);

    paper_save(h, pn.figures, maps(m).name);
    close(h);
end
