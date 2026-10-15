% PAPER_SUPP_FIG2_DENSITY_THICKNESS
%
% Supplementary Fig. 2: skull temperature rise vs skull thickness for the
% bottom ('low') and top ('high') skull density, from the density x
% thickness design of PRESTUS_thickness_04_density_sim.m (as plotted in
% PRESTUS_thickness_05_density_plot.m). The 2.5 mm benchmark skull
% (maximal transmission, see paper_fig4_thickness.m) is highlighted.
%
% Output: figures/paper/supp_fig2/suppfig2_density_thickness_riseT_skull.{pdf,png}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn.results = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a_density_thickness');
pn.figures = fullfile(rootpath, 'figures', 'paper', 'supp_fig2');

outcome      = 'riseT_skull';
highlight_mm = 2.5;

%% load (cf. PRESTUS_thickness_05_density_plot.m)

design = readtable(fullfile(pn.results, 'density_thickness_design.csv'), 'Delimiter', ',', 'TextType', 'string');
design.(outcome) = NaN(height(design), 1);
for i = 1:height(design)
    f = fullfile(pn.results, sprintf('sub-%03d', design.subject_id(i)), ...
        sprintf('sub-%03d_phantom%s.csv', design.subject_id(i), design.output_affix(i)));
    if ~isfile(f); continue; end
    R = readtable(f, 'Delimiter', ',');
    if ismember(outcome, R.Properties.VariableNames)
        design.(outcome)(i) = R.(outcome)(1);
    end
end
fprintf('%d/%d evaluations available.\n', sum(~isnan(design.(outcome))), height(design));

%% plot

S = paper_style();
levels = {'bottom', 'top'}; names = {'low', 'high'}; clrs = {S.low_line, S.high_line};

h = paper_figure([10.8 5.6]);
ax = axes(h, 'Units', 'centimeters', 'Position', [1.4 1.2 6.6 4.0]);
hold(ax, 'on');
y_max = 0.1 * ceil(10 * max(design.(outcome)) * 1.05);
patch(ax, highlight_mm + 0.35 * [-1 1 1 -1], [0 0 y_max y_max], [0.92 0.92 0.92], ...
    'EdgeColor', [0.4 0.4 0.4], 'LineStyle', '--');
for l = 1:numel(levels)
    sub = sortrows(design(design.skull_density_level == levels{l} & ~isnan(design.(outcome)), :), 'thickness_mm');
    plot(ax, sub.thickness_mm, sub.(outcome), '-o', 'Color', clrs{l}, 'LineWidth', S.linewidth, ...
        'MarkerSize', 3, 'MarkerFaceColor', clrs{l});
    [~, i_end] = max(sub.thickness_mm);
    text(ax, sub.thickness_mm(i_end) + 0.3, sub.(outcome)(i_end), ...
        sprintf('\\bf%s\\rm skull density', names{l}), 'Color', clrs{l});
end
set(ax, 'XLim', [min(design.thickness_mm) max(design.thickness_mm)], 'YLim', [0 y_max]);
xlabel(ax, 'Skull thickness (mm)');
ylabel(ax, 'Skull temperature rise (°C)');

paper_save(h, pn.figures, 'suppfig2_density_thickness_riseT_skull');
close(h);
