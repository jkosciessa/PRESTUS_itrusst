% PAPER_FIG4_THICKNESS
%
% Paper Fig. 4 (left): maximum brain intensity (Isppa_brain) as a function of skull
% thickness (thickness phantoms sub-301..312 / sub-351..361, reference
% medium, data/tussim/itrusst_protocol), marking the 2.5 mm skull used for
% the paramvar benchmark (sub-352; red line), the most transmissive
% thickness in the data (red dot) and the 6.5 mm default skull (grey, dashed).
%
% The transmission peak is a skull-layer resonance (thickness ~ lambda/2 in
% bone, 2.55 mm at the mid sound speed of 2550 m/s) and falls between the
% 0.5 mm thickness steps; with the current tissue parameters it lies at
% 2.5 mm, which is the paramvar benchmark (3 mm, used before, is ~4 % lower).
%
% Output: figures/paper/fig4/fig4a_isppa_brain_vs_thickness.{pdf,png}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn.data    = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol');
pn.figures = fullfile(rootpath, 'figures', 'paper', 'fig4');

default_mm   = 6.5;
benchmark_mm = 2.5; % thickness of the paramvar benchmark sub-352

%% load

[subjects, thickness] = itrusst_thickness_subjects(1:0.5:12);
Isppa_brain = nan(size(subjects));
for i = 1:numel(subjects)
    f = fullfile(pn.data, sprintf('sub-%03d', subjects(i)), sprintf('sub-%03d_phantom.csv', subjects(i)));
    if isfile(f)
        T = readtable(f);
        Isppa_brain(i) = T.Isppa_brain(1);
    else
        warning('File %s not found.', f);
    end
end

[I_max, i_max] = max(Isppa_brain);
I_bench = Isppa_brain(thickness == benchmark_mm);
fprintf('Most transmissive thickness: %.1f mm (%.2f W/cm^2); %.1f mm benchmark: %.2f W/cm^2 (%+.1f %%).\n', ...
    thickness(i_max), I_max, benchmark_mm, I_bench, 100 * (I_bench - I_max) / I_max);

%% plot

S = paper_style();
h = paper_figure(S.range_fig + [0.6 0]);
ax = axes(h, 'Units', 'centimeters', 'Position', S.range_axes + [0 0 0.6 0]);
hold(ax, 'on');
y_max = 5 * ceil(max(Isppa_brain) / 5);
plot(ax, [benchmark_mm benchmark_mm], [0 y_max], '-', 'Color', S.highlight, 'LineWidth', 1.5);
plot(ax, [default_mm default_mm], [0 y_max], '--', 'Color', S.ref_alt, 'LineWidth', 1.5);
plot(ax, thickness, Isppa_brain, 'k-', 'LineWidth', 1.5);
scatter(ax, thickness(i_max), I_max, 18, S.highlight, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
text(ax, benchmark_mm, y_max, sprintf('%g mm', benchmark_mm), 'Color', S.highlight, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax, default_mm, y_max, 'default', 'Color', S.ref_alt, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
set(ax, 'XLim', [min(thickness) max(thickness)], 'YLim', [0 y_max], 'XGrid', 'on', ...
    'GridColor', [0.85 0.85 0.85], 'GridAlpha', 1, 'Layer', 'bottom');
xlabel(ax, 'Skull thickness (mm)');
ylabel(ax, 'Brain I_{SPPA} (W/cm^2)');

paper_save(h, pn.figures, 'fig4a_isppa_brain_vs_thickness');
close(h);
