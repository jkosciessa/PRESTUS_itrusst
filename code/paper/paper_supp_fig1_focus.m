% PAPER_SUPP_FIG1_FOCUS
%
% Supplementary Fig. 1: focal position and volume under acoustic property
% variation (sub-002, skull benchmark, itrusst_protocol1a).
%   (a) range panels of the focal distance and the -3 dB (FWHM) focal
%       volume for the acoustic properties of brain and skull
%   (b) on-axis I_SPPA for the bottom ('low') and top ('high') skull sound
%       speed, with the skull boundaries dashed, and the -3 dB focal masks
%       of both runs (low / high / both)
%
% Output: figures/paper/supp_fig1/*.{pdf,png}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn.data    = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a', 'sub-002');
pn.simnibs = fullfile(rootpath, 'data', 'simnibs');
pn.figures = fullfile(rootpath, 'figures', 'paper', 'supp_fig1');

S = paper_style();

%% (a) focal distance and volume

symbols = {'c', '\rho', '\alpha'};
[b, ~] = load_paramvar_results(pn.data, 2, {'brain'}, {'_b'}, ...
    {'sound_speed', 'density', 'alpha_coeff'}, {'_c_', '_rho_', '_a0_'}, 1:2, symbols);
[s, ~] = load_paramvar_results(pn.data, 2, {'skull'}, {'_s'}, ...
    {'sound_speed', 'density', 'alpha_coeff'}, {'_c_', '_rho_', '_a0_'}, 1:2, symbols);
tissues = [repmat({'b'}, 1, 3), repmat({'s'}, 1, 3)];

metrics = {'real_focal_distance_mm', 'Focal distance (mm)',  'suppfig1a_focal_distance'; ...
           'focal_vol_3dB_mm3',      'Focal volume (mm^3)',  'suppfig1a_focal_volume'};
for i = 1:size(metrics, 1)
    data = [b.(metrics{i, 1}), s.(metrics{i, 1})];
    lims = [floor(min(data(:))) - 1, ceil(max(data(:))) + 1];
    h = paper_figure(S.range_fig);
    ax = axes(h, 'Units', 'centimeters', 'Position', S.range_axes);
    paper_range_panel(ax, data, [symbols symbols], tissues, metrics{i, 2}, lims);
    paper_save(h, pn.figures, metrics{i, 3});
    close(h);
end
vol = [b.focal_vol_3dB_mm3, s.focal_vol_3dB_mm3];

%% (b) skull sound speed: axial profiles and focal masks

tissue_info = niftiinfo(fullfile(pn.simnibs, 'm2m_sub-002', 'final_tissues.nii.gz'));
tissue = round(double(niftiread(tissue_info)));
res = tissue_info.PixelDimensions(1);
cx = round(size(tissue, 1) / 2);
z_mm = ((1:size(tissue, 2)) - 1) * res;
x_mm = ((1:size(tissue, 1)) - cx) * res;
skull_z = z_mm(tissue(cx, :) == 4);

h = paper_figure([13.0 5.4]);
ax = axes(h, 'Units', 'centimeters', 'Position', [1.2 1.1 6.8 3.8]);
hold(ax, 'on');
names = {'low', 'high'}; clrs = {S.low_line, S.high_line};
y_max = 0;
masks = cell(1, 2);
for l = 1:2
    I = double(niftiread(fullfile(pn.data, 'nii', sprintf('sub-002_phantom_T1w_s_c_%d_Ipa.nii.gz', l))));
    plot(ax, z_mm, I(cx, :), '-', 'Color', clrs{l}, 'LineWidth', S.linewidth);
    y_max = max(y_max, max(I(cx, :)));
    masks{l} = double(niftiread(fullfile(pn.data, 'nii', sprintf('sub-002_phantom_T1w_s_c_%d_focal_mask_3dB.nii.gz', l)))) > 0;
end
y_max = 2 * ceil(1.2 * y_max / 2);
for z = [min(skull_z) - res/2, max(skull_z) + res/2]
    plot(ax, [z z], [0 y_max], '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1);
end
set(ax, 'XLim', [0 max(z_mm)], 'YLim', [0 y_max]);
xlabel(ax, 'Distance from transducer (mm)');
ylabel(ax, 'Intensity (W/cm^2)');
text(ax, 0.58 * max(z_mm), 0.95 * y_max, '\bflow\rm skull sound speed', 'Color', S.low_line, 'FontSize', S.fontsize - 1);
text(ax, 0.58 * max(z_mm), 0.86 * y_max, '\bfhigh\rm skull sound speed', 'Color', S.high_line, 'FontSize', S.fontsize - 1);

% focal masks: low only / high only / both
axm = axes(h, 'Units', 'centimeters', 'Position', [10.2 1.1 2.4 3.8]);
category = zeros(size(masks{1}));
category(masks{1} & ~masks{2}) = 1;
category(~masks{1} & masks{2}) = 2;
category(masks{1} & masks{2}) = 3;
imagesc(axm, x_mm, z_mm, category', [0 3]);
colormap(axm, [0.93 0.93 0.93; S.low_line; S.high_line; 0.55 0.35 0.75]);
set(axm, 'YDir', 'reverse', 'XLim', [-30 30], 'YLim', [0 max(z_mm)], 'XTick', [-30 0 30]);
xlabel(axm, 'Lateral (mm)');
ylabel(axm, 'Distance from transducer (mm)');
text(axm, -27, 6, sprintf('-3 dB focal volume\n%.0f-%.0f mm^3', min(vol(:)), max(vol(:))), ...
    'FontAngle', 'italic', 'FontSize', S.fontsize - 1, 'VerticalAlignment', 'top');

paper_save(h, pn.figures, 'suppfig1b_skull_soundspeed_profiles');
close(h);
