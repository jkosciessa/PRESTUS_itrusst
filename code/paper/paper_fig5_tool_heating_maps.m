% PAPER_FIG5_TOOL_HEATING_MAPS
%
% Paper Fig. 5 (candidate): heat-rise maps of the emulated toolbox defaults
% (PRESTUS, k-Plan*, BabelBrain*; data/tussim/medium_<tool>) in the plane
% through the beam axis of the 3D benchmarks (sub-102: skull, sub-101: no
% skull), with the skull outline.
%   top:    maximum temperature rise per tool (shared colour range)
%   bottom: pairwise differences (k-Plan* - PRESTUS, BabelBrain* - PRESTUS,
%           k-Plan* - BabelBrain*), colour range symmetric around zero and
%           shared across pairs (black: cooler, red: warmer), as in
%           paper_fig3_heating_difference_maps.m.
% All tools are simulated with PRESTUS/k-Wave (axisymmetric acoustics,
% thermal on the radially expanded 3D field); only the medium differs. The
% heat-rise NIfTIs hold the plane through the beam axis.
%
% Output: figures/paper/fig5/fig5_heatmaps_tools_{skull,noskull}.{pdf,png}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn.tussim  = fullfile(rootpath, 'data', 'tussim');
pn.simnibs = fullfile(rootpath, 'data', 'simnibs');
pn.figures = fullfile(rootpath, 'figures', 'paper', 'fig5');

%% settings

tools  = {'PRESTUS', 'kPlan', 'BabelBrain'};
labels = {'PRESTUS', 'k-Plan*', 'BabelBrain*'};
pairs  = [2 1; 3 1; 2 3]; % minuend, subtrahend (indices into tools)
benchmarks = struct( ...
    'sub',  {'sub-102', 'sub-101'}, ...
    'name', {'skull',   'noskull'});
depth_lim   = [0 90];   % mm from the transducer
lateral_lim = [-35 35]; % mm from the axis
axis_vox    = 70;       % beam axis (lateral voxel)

S = paper_style();
cmap_diff = [linspace(0, 1, 128)', linspace(0, 1, 128)', linspace(0, 1, 128)'; ... % black -> white
             ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)'];            % white -> red
cmap_rise = cmap_diff(129:end, :);                                                 % white -> red

%% plot

for b = 1:numel(benchmarks)
    sub = benchmarks(b).sub;

    tissue = round(niftiread(fullfile(pn.simnibs, ['m2m_' sub], 'final_tissues.nii.gz')));
    skull = double(squeeze(tissue(:, axis_vox, :)) == 4);

    rise = cell(1, numel(tools));
    for t = 1:numel(tools)
        f = fullfile(pn.tussim, ['medium_' tools{t}], sub, 'nii', [sub '_phantom_T1w_heatrise.nii.gz']);
        rise{t} = double(niftiread(f)); % plane through the beam axis [lateral x axial]
        fprintf('%s %s: max rise %.3f °C\n', sub, tools{t}, max(rise{t}(:)));
    end
    info = niftiinfo(f);
    res = info.PixelDimensions(1);
    x_mm = ((1:size(rise{1}, 1)) - axis_vox) * res;
    z_mm = ((1:size(rise{1}, 2)) - 1) * res;

    rise_lim = ceil(max(cellfun(@(r) max(r(:)), rise)) * 10) / 10;
    diffs = arrayfun(@(p) rise{pairs(p, 1)} - rise{pairs(p, 2)}, 1:size(pairs, 1), 'UniformOutput', false);
    diff_lim = ceil(max(cellfun(@(d) max(abs(d(:))), diffs)) * 10) / 10;
    for p = 1:size(pairs, 1)
        fprintf('%s %s - %s in [%.3f, %.3f] °C\n', sub, tools{pairs(p, 1)}, tools{pairs(p, 2)}, ...
            min(diffs{p}(:)), max(diffs{p}(:)));
    end

    ax_w = 3.6;
    ax_h = ax_w * diff(depth_lim) / diff(lateral_lim);
    h = paper_figure([3 * (ax_w + 0.3) + 1.6, 2 * (ax_h + 0.7) + 0.2]);
    for row = 1:2
        for col = 1:3
            pos = [0.3 + (col - 1) * (ax_w + 0.3), 0.3 + (2 - row) * (ax_h + 0.7), ax_w, ax_h];
            ax = axes(h, 'Units', 'centimeters', 'Position', pos);
            if row == 1
                imagesc(ax, x_mm, z_mm, rise{col}', [0 rise_lim]);
                colormap(ax, cmap_rise);
                ttl = labels{col};
            else
                imagesc(ax, x_mm, z_mm, diffs{col}', [-diff_lim diff_lim]);
                colormap(ax, cmap_diff);
                ttl = sprintf('%s - %s', labels{pairs(col, 1)}, labels{pairs(col, 2)});
            end
            hold(ax, 'on');
            if any(skull(:))
                contour(ax, x_mm, z_mm, skull', [0.5 0.5], 'Color', [0.2 0.2 0.2], 'LineWidth', 0.5);
            end
            axis(ax, 'image');
            set(ax, 'XLim', lateral_lim, 'YLim', depth_lim, 'YDir', 'reverse', 'Visible', 'off');
            title(ax, ttl, 'Visible', 'on');
            if col == 3
                cb = colorbar(ax, 'Units', 'centimeters', ...
                    'Position', [pos(1) + ax_w + 0.3, pos(2) + 0.2 * ax_h, 0.25, 0.6 * ax_h]);
                if row == 1
                    cb.Ticks = [0 rise_lim];
                    cb.Label.String = 'Heat rise (°C)';
                else
                    cb.Ticks = [-diff_lim 0 diff_lim];
                    cb.Label.String = 'Heating difference (°C)';
                end
                cb.TickDirection = 'out';
                cb.LineWidth = 0.5;
                set(ax, 'Position', pos);
            end
        end
    end

    paper_save(h, pn.figures, ['fig5_heatmaps_tools_' benchmarks(b).name]);
    close(h);
end
