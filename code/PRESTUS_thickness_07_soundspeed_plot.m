% PRESTUS_THICKNESS_07_SOUNDSPEED_PLOT
%
% Plots the skull sound speed x thickness design of
% PRESTUS_thickness_06_soundspeed_sim.m: each outcome (Isppa_brain,
% Ipa_target, riseT_brain, riseT_skull) vs. c_skull, one panel per skull
% thickness, with the literature bottom / mid / top marked.
%
% Brain I_SPPA panels overlay a plane-wave layer estimate (water | skull |
% brain at normal incidence, mid density and attenuation), scaled to the
% simulation at mid c_skull: it predicts where the transmission dips and
% peaks fall, not their absolute size (no focusing, no oblique incidence).
% Missing or incomplete result files are reported and left as gaps.
%
% Output: figures/thickness/soundspeed_thickness_<outcome>.{png,eps}

%% ── Paths ───────────────────────────────────────────────────────────────────

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'));
rootpath = pwd;

variant_name = 'itrusst_protocol1a_soundspeed_thickness';   % must match PRESTUS_thickness_06
pn.results = fullfile(rootpath, 'data', 'tussim', variant_name);
pn.figures = fullfile(rootpath, 'figures', 'thickness');
if ~isfolder(pn.figures); mkdir(pn.figures); end

outcomes = {'Isppa_brain', 'Ipa_target', 'riseT_brain', 'riseT_skull'};
outcome_labels = {'Brain I_{SPPA} (W/cm^2)', 'I_{PA} at target (W/cm^2)', ...
                  'Brain temperature rise (°C)', 'Skull temperature rise (°C)'};
freq_hz = 500e3;   % transducer frequency (config_itrusst_protocol.yaml)

%% ── Load the design table and join results ─────────────────────────────────

design = readtable(fullfile(pn.results, 'soundspeed_thickness_design.csv'), 'Delimiter', ',', 'TextType', 'string');
n_rows = height(design);

for o = 1:numel(outcomes)
    design.(outcomes{o}) = NaN(n_rows, 1);
end
design.available = false(n_rows, 1);

n_missing = 0;
for i = 1:n_rows
    result_file = fullfile(pn.results, sprintf('sub-%03d', design.subject_id(i)), ...
        sprintf('sub-%03d_phantom%s.csv', design.subject_id(i), design.output_affix(i)));
    if ~isfile(result_file)
        n_missing = n_missing + 1;
        continue
    end
    R = readtable(result_file, 'Delimiter', ',');
    if ~all(ismember(outcomes, R.Properties.VariableNames))
        n_missing = n_missing + 1;   % ran but incomplete (e.g. killed mid-thermal-stage)
        continue
    end
    design.available(i) = true;
    for o = 1:numel(outcomes)
        design.(outcomes{o})(i) = R.(outcomes{o})(1);
    end
end

fprintf('Sound speed x thickness: %d/%d evaluations available (%d missing/incomplete).\n', ...
    n_rows - n_missing, n_rows, n_missing);

%% ── Literature bounds + plane-wave layer estimate ──────────────────────────

P = readtable(fullfile(rootpath, 'data', 'tissue_parameters.csv'));
prop = @(tissue, level, name) P.(name)(strcmp(P.tissue, tissue) & strcmp(P.level, level));
c_bounds = [prop('skull', 'bottom', 'sound_speed'), prop('skull', 'mid', 'sound_speed'), ...
            prop('skull', 'top', 'sound_speed')];

layer.Z_water = prop('water', 'mid', 'density') * prop('water', 'mid', 'sound_speed');
layer.Z_brain = prop('brain', 'mid', 'density') * prop('brain', 'mid', 'sound_speed');
layer.rho     = prop('skull', 'mid', 'density');
% dB/(MHz^y cm) -> Np/m at freq_hz
layer.alpha   = prop('skull', 'mid', 'alpha_coeff') * (freq_hz / 1e6)^prop('skull', 'mid', 'alpha_power') ...
                * 100 / (20 * log10(exp(1)));
layer.freq    = freq_hz;

%% ── Plot: one figure per outcome, one panel per thickness ──────────────────

thicknesses = unique(design.thickness_mm);
c_fine = linspace(min(design.skull_sound_speed), max(design.skull_sound_speed), 200);

for o = 1:numel(outcomes)
    h = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 6 * numel(thicknesses) 7], ...
        'InvertHardcopy', 'off', 'Renderer', 'painters');
    tl = tiledlayout(h, 1, numel(thicknesses), 'TileSpacing', 'compact');
    for t = 1:numel(thicknesses)
        ax = nexttile(tl);
        hold(ax, 'on');
        sub = sortrows(design(design.thickness_mm == thicknesses(t) & design.available, :), ...
            'skull_sound_speed');
        xline(ax, c_bounds, ':', {'bottom', 'mid', 'top'}, 'Color', [0.5 0.5 0.5], ...
            'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
        if strcmp(outcomes{o}, 'Isppa_brain') && ~isempty(sub)
            i_mid = find(sub.skull_sound_speed == c_bounds(2), 1);
            if ~isempty(i_mid)
                T_fine = layer_transmission(c_fine, thicknesses(t) * 1e-3, layer);
                T_mid  = layer_transmission(c_bounds(2), thicknesses(t) * 1e-3, layer);
                plot(ax, c_fine, T_fine / T_mid * sub.Isppa_brain(i_mid), '-', ...
                    'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            end
        end
        plot(ax, sub.skull_sound_speed, sub.(outcomes{o}), '-o', 'Color', [0.7 0 0], ...
            'MarkerFaceColor', [0.7 0 0], 'MarkerSize', 4, 'LineWidth', 1.2);
        title(ax, sprintf('%.1f mm skull', thicknesses(t)), 'FontWeight', 'normal');
        xlabel(ax, 'c_{skull} (m/s)');
        if t == 1; ylabel(ax, outcome_labels{o}, 'Interpreter', 'tex'); end
        set(ax, 'Box', 'off', 'FontSize', 9, 'XLim', [min(c_fine) - 30, max(c_fine) + 30]);
    end
    if strcmp(outcomes{o}, 'Isppa_brain')
        title(tl, 'red: simulation; grey: plane-wave layer estimate, scaled at mid', 'FontSize', 9);
    end
    fig_name = sprintf('soundspeed_thickness_%s', outcomes{o});
    saveas(h, fullfile(pn.figures, fig_name), 'epsc');
    saveas(h, fullfile(pn.figures, fig_name), 'png');
    fprintf('  saved %s\n', fullfile(pn.figures, fig_name));
end

%% ── Local functions ──────────────────────────────────────────────────────

function T = layer_transmission(c, L, layer)
% Intensity transmission water -> skull layer (thickness L, sound speed c)
% -> brain, plane wave at normal incidence (transfer matrix).
T = nan(size(c));
for i = 1:numel(c)
    k = 2 * pi * layer.freq / c(i) - 1i * layer.alpha;
    Z = layer.rho * c(i);
    M = [cos(k * L), 1i * Z * sin(k * L); 1i * sin(k * L) / Z, cos(k * L)];
    t_p = 2 / (M(1, 1) + M(1, 2) / layer.Z_brain + layer.Z_water * M(2, 1) + ...
               layer.Z_water * M(2, 2) / layer.Z_brain);
    T(i) = abs(t_p)^2 * layer.Z_water / layer.Z_brain;
end
end
