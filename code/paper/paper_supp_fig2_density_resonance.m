% PAPER_SUPP_FIG2_DENSITY_RESONANCE
%
% Supplementary Fig. 2 (schematic): why a denser skull heats less at most
% thicknesses but more at the half-wave resonance (2.5 mm), as seen in
% paper_supp_fig2_density_thickness.m.
%
% 1D plane-wave model at normal incidence: water | skull plate | brain,
% 500 kHz, skull sound speed at its mid level (2550 m/s, lambda/2 = 2.55 mm),
% skull density at its bottom/top level (1080 / 1991 kg/m^3), acoustic
% attenuation and heat absorption at their mid levels; all values from
% data/tissue_parameters.csv as used in PRESTUS_thickness_04_density_sim.m.
% Heat deposition as in PRESTUS (thermal_simulation.m): Q = alpha * p^2 / (rho * c).
%
% Mechanism: at d = n * lambda/2 the plate is transparent whatever its
% impedance, but the standing wave inside it grows with the impedance
% contrast m = Z_skull / Z_water. Without loss, the mean heat deposition in
% the plate scales with (m + 1/m), so a denser skull absorbs more. Away from
% resonance the denser plate reflects more and absorbs less.
%
% (a) schematic: pressure amplitude |p(x)| across the layers for the two
%     benchmark skulls, 2.5 mm (~lambda/2, resonant) and 6.5 mm (~5/4 lambda,
%     anti-resonant), with the high/low ratio of the peak heat deposition
% (b) model prediction: high/low ratio of the peak heat deposition in the
%     skull vs thickness (1D model; axis symmetric around 1)
% (c) heating data: skull temperature rise vs thickness for low/high density
%     (PRESTUS riseT_skull, as in paper_supp_fig2_density_thickness.m). The
%     curves cross only at the first resonance: the swept density also sets
%     the volumetric heat capacity rho * c_p in kWaveDiffusion, so a denser
%     skull warms less per unit heat.
%
% Limits: 1D, normal incidence, no thermal model. The simulated focused
% beam meets a curved shell, so its resonances fall at slightly larger
% thicknesses than n * lambda/2 and persist in thick skulls.
%
% Output: figures/paper/supp_fig2/suppfig2_density_combined.{pdf,png}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn.results = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a_density_thickness');
pn.figures = fullfile(rootpath, 'figures', 'paper', 'supp_fig2');

%% medium (data/tissue_parameters.csv)

P = readtable(fullfile(rootpath, 'data', 'tissue_parameters.csv'), 'TextType', 'string');
prop = @(t, l, f) P.(f)(P.tissue == t & P.level == l);

f0 = 500e3;                                  % [Hz], config_itrusst_protocol.yaml
w  = 2 * pi * f0;
db2np = 1 / (20 * log10(exp(1)));            % dB -> Np

c1 = prop("water", "mid", "sound_speed"); rho1 = prop("water", "mid", "density");
c3 = prop("brain", "mid", "sound_speed"); rho3 = prop("brain", "mid", "density");
c2 = prop("skull", "mid", "sound_speed");
rho2 = [prop("skull", "bottom", "density"), prop("skull", "top", "density")];
alpha_ac = prop("skull", "mid", "alpha_coeff")      * (f0 / 1e6) * 100 * db2np; % [Np/m], alpha_power = 1
alpha_h  = prop("skull", "mid", "absorption_coeff") * (f0 / 1e6) * 100 * db2np; % [Np/m], heat

lambda2 = c2 / f0 * 1e3;                     % wavelength in skull [mm]
Z1 = rho1 * c1; Z3 = rho3 * c3;
fprintf('lambda/2 in skull: %.2f mm; m = Z_skull/Z_water: %.2f (low), %.2f (high)\n', ...
    lambda2 / 2, rho2 * c2 / Z1);

%% 1D layer model

thick_mm = 0.5:0.01:12;
q_peak = nan(numel(thick_mm), 2);   % maximum heat deposition within the skull
for l = 1:2
    for i = 1:numel(thick_mm)
        [~, q_peak(i, l)] = layer_heat(thick_mm(i) * 1e-3, w, c2, rho2(l), alpha_ac, alpha_h, Z1, Z3);
    end
end

%% PRESTUS density x thickness design (cf. paper_supp_fig2_density_thickness.m)

outcomes = {'riseT_skull'};
design = readtable(fullfile(pn.results, 'density_thickness_design.csv'), 'Delimiter', ',', 'TextType', 'string');
for o = 1:numel(outcomes); design.(outcomes{o}) = NaN(height(design), 1); end
for i = 1:height(design)
    f = fullfile(pn.results, sprintf('sub-%03d', design.subject_id(i)), ...
        sprintf('sub-%03d_phantom%s.csv', design.subject_id(i), design.output_affix(i)));
    if ~isfile(f); continue; end
    R = readtable(f, 'Delimiter', ',');
    for o = 1:numel(outcomes)
        if ismember(outcomes{o}, R.Properties.VariableNames)
            design.(outcomes{o})(i) = R.(outcomes{o})(1);
        end
    end
end
lo = sortrows(design(design.skull_density_level == "bottom", :), 'thickness_mm');
hi = sortrows(design(design.skull_density_level == "top", :), 'thickness_mm');

%% plot

S = paper_style();
M = struct('w', w, 'c1', c1, 'c2', c2, 'rho2', rho2, 'alpha_ac', alpha_ac, ...
    'alpha_h', alpha_h, 'Z1', Z1, 'Z3', Z3, 'lambda2', lambda2);
r_model = q_peak(:, 2) ./ q_peak(:, 1);

% Supp. Fig. 2: (a) schematic, (b) model prediction, (c) heating data
h = paper_figure([17.4 11.0]);
draw_profiles(h, [1.9 7.0 6.4 3.0; 10.4 7.0 6.4 3.0], M, S);
draw_model(h, [1.9 1.2 6.4 4.0], thick_mm, r_model, lambda2, S);
draw_data(h, [10.4 1.2 6.4 4.0], lo, hi, lambda2, S);
ann = axes(h, 'Units', 'centimeters', 'Position', [0 0 h.Position(3:4)], 'Visible', 'off', ...
    'XLim', [0 h.Position(3)], 'YLim', [0 h.Position(4)]);
text(ann, 0.1, 10.9, 'a', 'FontWeight', 'bold', 'FontSize', 12, 'VerticalAlignment', 'top');
text(ann, 0.1, 5.8, 'b', 'FontWeight', 'bold', 'FontSize', 12, 'VerticalAlignment', 'top');
text(ann, 8.9, 5.8, 'c', 'FontWeight', 'bold', 'FontSize', 12, 'VerticalAlignment', 'top');
paper_save(h, pn.figures, 'suppfig2_density_combined');
close(h);

%% helpers

function draw_model(h, pos, thick_mm, r_model, lambda2, S)
% 1D model: high/low ratio of the peak heat deposition in the skull
ax = axes(h, 'Units', 'centimeters', 'Position', pos);
hold(ax, 'on');
r_lim = [0.5 1.5];   % symmetric around a ratio of 1
patch(ax, [1 12 12 1], [1 1 r_lim(2) r_lim(2)], [0.95 0.95 0.95], 'EdgeColor', 'none');
text(ax, 11.8, r_lim(2) - 0.02, 'denser skull absorbs more', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'Color', S.ref_alt);
text(ax, 11.8, r_lim(1) + 0.02, 'denser skull absorbs less', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', 'Color', S.ref_alt);
plot_resonances(ax, lambda2, r_lim, S);
plot(ax, [1 12], [1 1], '-', 'Color', S.ref, 'LineWidth', 0.5);
plot(ax, thick_mm, r_model, '-', 'Color', S.ref, 'LineWidth', S.linewidth);
set(ax, 'XLim', [1 12], 'YLim', r_lim, 'YTick', 0.5:0.25:1.5, 'Layer', 'top');
xlabel(ax, 'Skull thickness (mm)');
ylabel(ax, {'Peak heat in skull (1D model)', 'high / low density'});
end

function draw_data(h, pos, lo, hi, lambda2, S)
% PRESTUS skull temperature rise vs thickness for low/high density (as in
% paper_supp_fig2_density_thickness.m), 2.5 mm benchmark skull highlighted
ax = axes(h, 'Units', 'centimeters', 'Position', pos);
hold(ax, 'on');
sets = {lo, hi}; clrs = {S.low_line, S.high_line};
y_max = 0.4;
patch(ax, 2.5 + 0.35 * [-1 1 1 -1], [0 0 y_max y_max], [0.92 0.92 0.92], ...
    'EdgeColor', [0.4 0.4 0.4], 'LineStyle', '--');
plot_resonances(ax, lambda2, [0 y_max], S);
for l = 1:2
    sub = sets{l}(~isnan(sets{l}.riseT_skull), :);
    plot(ax, sub.thickness_mm, sub.riseT_skull, '-o', 'Color', clrs{l}, 'LineWidth', S.linewidth, ...
        'MarkerSize', 3, 'MarkerFaceColor', clrs{l});
end
text(ax, 11.8, 0.395, '\bflow\rm skull density', 'Color', clrs{1}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
text(ax, 11.8, 0.365, '\bfhigh\rm skull density', 'Color', clrs{2}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
set(ax, 'XLim', [1 12], 'YLim', [0 y_max], 'Layer', 'top');
xlabel(ax, 'Skull thickness (mm)');
ylabel(ax, 'Skull temperature rise (°C)');
end

function draw_profiles(h, pos, M, S)
% |p(x)| at a resonant and an anti-resonant thickness; pos = [2 x 4] axes
% positions (cm), one row per case
clrs = {S.low_line, S.high_line};
skull_patch = [0.96 0.88 0.88];
cases = [2.5 6.5];                           % [mm], the benchmark skulls: ~lambda/2 and ~5/4 lambda
case_names = {sprintf('2.5 mm skull: resonant (\\approx\\lambda/2 = %.2f mm)', M.lambda2 / 2), ...
              sprintf('6.5 mm skull: anti-resonant (\\approx5\\lambda/4 = %.2f mm)', 5 * M.lambda2 / 4)};
case_notes = {'denser skull: stronger standing wave', 'denser skull: more reflected'};
for k = 1:2
    d = cases(k);
    x = linspace(-5, d + 5, 800);            % [mm], skull from 0 to d
    ax = axes(h, 'Units', 'centimeters', 'Position', pos(k, :));
    hold(ax, 'on');
    patch(ax, [0 d d 0], [0 0 4 4], skull_patch, 'EdgeColor', 'none');
    q = zeros(1, 2);
    for l = 1:2
        [~, q(l), px] = layer_heat(d * 1e-3, M.w, M.c2, M.rho2(l), M.alpha_ac, M.alpha_h, M.Z1, M.Z3, x * 1e-3, M.c1);
        plot(ax, x, abs(px), '-', 'Color', clrs{l}, 'LineWidth', S.linewidth);
    end
    plot(ax, x([1 end]), [1 1], ':', 'Color', S.ref_alt, 'LineWidth', 0.5);
    set(ax, 'XLim', x([1 end]), 'YLim', [0 4], 'YTick', 0:1:4, 'XTick', -4:2:12);
    text(ax, 0, 1.02, case_names{k}, 'Units', 'normalized', 'VerticalAlignment', 'bottom');
    text(ax, x(end) - 0.2, 3.5, {case_notes{k}, sprintf('heat in skull, high/low: %.2f', q(2) / q(1))}, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    if k == 1
        text(ax, x(end) - 0.2, 2.1, '\bflow\rm density', 'Color', clrs{1}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        text(ax, x(end) - 0.2, 1.65, '\bfhigh\rm density', 'Color', clrs{2}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        text(ax, -2.5, 0.15, 'water', 'Color', S.ref_alt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        text(ax, d / 2, 0.15, 'skull', 'Color', S.skull, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        text(ax, d + 2.5, 0.15, 'brain', 'Color', S.ref_alt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    xlabel(ax, 'Depth from outer skull surface (mm)');
    ylabel(ax, {'Pressure amplitude', '(incoming wave = 1)'});
end
end

function plot_resonances(ax, lambda2, y_lim, S)
% dashed lines at the half-wave resonances n * lambda/2
for n = 1:floor(12 / (lambda2 / 2))
    x_n = n * lambda2 / 2;
    plot(ax, [x_n x_n], y_lim, '--', 'Color', S.ref_alt, 'LineWidth', 0.5);
    text(ax, x_n, y_lim(2), sprintf('%d\\lambda/2', n), 'Color', S.ref_alt, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
end

function [frac, q_max, px] = layer_heat(d, w, c2, rho2, alpha_ac, alpha_h, Z1, Z3, x, c1)
% Plane wave (unit amplitude, normal incidence) on a lossy layer of
% thickness d between half-spaces of impedance Z1 (front) and Z3 (back).
% frac: power absorbed in the layer per unit area, Q = alpha_h |p|^2/(rho c),
% relative to the incident intensity 1/(2 Z1); q_max: its maximum within the
% layer, scaled alike [1/m]. px: |p| at positions x [m]
% (layer from 0 to d), c1 = sound speed in front. Convention exp(j(wt - kx)).
k2 = w / c2 - 1i * alpha_ac;
Z2 = w * rho2 / k2;
r23 = (Z3 - Z2) / (Z3 + Z2);
% inside: p = A exp(-j k2 x) + B exp(j k2 x), B = r23 A exp(-2j k2 d)
b = r23 * exp(-2i * k2 * d);
Zin = Z2 * (1 + b) / (1 - b);
R = (Zin - Z1) / (Zin + Z1);
A = (1 + R) / (1 + b);
B = b * A;
xi = linspace(0, d, 400);
p_in = A * exp(-1i * k2 * xi) + B * exp(1i * k2 * xi);
q_in = alpha_h * abs(p_in).^2 / (rho2 * c2) * 2 * Z1;
frac = trapz(xi, q_in);
q_max = max(q_in);
if nargin < 9; px = []; return; end
px = zeros(size(x));
pre = x < 0; ins = x >= 0 & x <= d; post = x > d;
% water: incident + reflected; brain: transmitted (lossless outside)
px(pre)  = abs(exp(-1i * w / c1 * x(pre)) + R * exp(1i * w / c1 * x(pre)));
px(ins)  = abs(A * exp(-1i * k2 * x(ins)) + B * exp(1i * k2 * x(ins)));
p_d = A * exp(-1i * k2 * d) + B * exp(1i * k2 * d);
px(post) = abs(p_d) * ones(1, sum(post));
end
