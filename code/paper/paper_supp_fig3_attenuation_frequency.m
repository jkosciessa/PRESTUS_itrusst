% PAPER_SUPP_FIG3_ATTENUATION_FREQUENCY
%
% Supplementary Fig. 3: reported skull attenuation as a function of
% frequency, in dB/cm. (a) whole skull or layer not separated, (b) cortical
% tables and diploë. Values are the measured (or study-fitted) values per
% frequency taken from each paper, not values rescaled with a power law; see
% data/literature_skull-attenuation_values.csv (source location, extraction
% and caveats per row; rows with plot = 0 are kept for reference only; the
% s1_row column links rows to Supplementary Table S1).
%
% Conversion to dB/cm: 1 Np = 8.686 dB. Insertion losses (dB, not
% reflection-corrected) are divided by the reported skull thickness.
% Values that include reflection (insertion loss per thickness, Lee 2020,
% and Riis 2022, whose reflection term is a fitted intercept) are flagged in
% the reflection_included column and drawn with open markers and dotted
% lines, as they are not attenuation coefficients. Where a paper reports one
% value per skull or location (Ammi 2008, Riis 2022, Albelda Gimeno 2019)
% the mean ± SD across skulls/locations is computed here. Series and legend
% entries follow the label column of the data file. Error bars elsewhere are
% as reported (SD, or ± of undefined type for White 2006 and Ammi 2008).
%
% Background: the shaded band marks frequencies common in neuromodulation
% (0.2-0.65 MHz); the dotted line marks 0.5 MHz, the simulation frequency.
% Horizontal bars at 0.5 MHz mark the skull attenuation LOW / reference /
% HIGH of Table 1, read from data/tissue_parameters.csv.
%
% Output: figures/paper/supp_fig3/suppfig3_skull_attenuation_frequency.{pdf,png}
%         (PDF is vector, see paper_save) and the plotted values (dB/cm)
%         as ..._plotted.csv

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn.data    = fullfile(rootpath, 'data', 'literature_skull-attenuation_values.csv');
pn.bounds  = fullfile(rootpath, 'data', 'tissue_parameters.csv');
pn.figures = fullfile(rootpath, 'figures', 'paper', 'supp_fig3');

np2db  = 20 / log(10);                                     % 8.686 dB per Np
pooled = {'ammi2008', 'riis2022', 'gimeno2019'};           % one value per skull/location
xl     = [0 2.3];                                          % MHz
fband  = [0.2 0.65];                                       % common neuromodulation frequencies (MHz)
fsim   = 0.5;                                              % simulation frequency (MHz)

%% load and convert to dB/cm

numvars = {'freq_mhz', 'value', 'sd', 'thickness_mm', 'thickness_sd_mm', 'plot', 'reflection_included'};
opts = detectImportOptions(pn.data, 'Delimiter', ',', 'TextType', 'char');
opts = setvartype(opts, numvars, 'double');
opts = setvartype(opts, setdiff(opts.VariableNames, numvars), 'char');
T = readtable(pn.data, opts);
T = T(T.plot == 1, :);

T.db_cm    = NaN(height(T), 1);
T.db_cm_sd = NaN(height(T), 1);
for i = 1:height(T)
    switch T.quantity{i}
        case {'attenuation_coefficient', 'insertion_loss_per_length'}
            switch T.unit{i}
                case 'Np/m';  k = np2db / 100;
                case 'Np/cm'; k = np2db;
                case 'dB/cm'; k = 1;
                otherwise; error('unit %s not handled (row %d)', T.unit{i}, i);
            end
        case 'insertion_loss'
            assert(strcmp(T.unit{i}, 'dB') && ~isnan(T.thickness_mm(i)), 'row %d: needs dB and a thickness', i);
            k = 10 / T.thickness_mm(i);                     % dB -> dB/cm
        otherwise
            error('quantity %s not handled (row %d)', T.quantity{i}, i);
    end
    T.db_cm(i)    = abs(T.value(i)) * k;
    T.db_cm_sd(i) = T.sd(i) * k;
end
T.il = T.reflection_included == 1;

%% Table 1 skull attenuation bounds at 0.5 MHz (dB/cm)

B = readtable(pn.bounds, 'TextType', 'char');
B = B(strcmp(B.tissue, 'skull'), :);
lvl = {'bottom', 'mid', 'top'};
bnd = NaN(1, 3);
for il = 1:3
    r = B(strcmp(B.level, lvl{il}), :);
    bnd(il) = r.alpha_coeff * fsim ^ r.alpha_power;        % PRESTUS stores the 1 MHz coefficient
end
fprintf('  Table 1 skull attenuation at %.1f MHz: LOW %.2f, reference %.2f, HIGH %.2f dB/cm\n', fsim, bnd);

%% series: one line per label (pooled studies: mean ± SD across skulls/locations)

T.series = T.label;

P = table();
for s = unique(T.series, 'stable')'
    Ts = T(strcmp(T.series, s{1}), :);
    for f = unique(Ts.freq_mhz)'
        Tf = Ts(Ts.freq_mhz == f, :);
        r = Tf(1, {'study', 'ref_key', 'specimen', 'layer', 'series', 'label', 'il'});
        r.freq_mhz = f;
        if height(Tf) > 1
            assert(ismember(Tf.ref_key{1}, pooled), 'several values for %s at %g MHz', s{1}, f);
            r.db_cm = mean(Tf.db_cm); r.db_cm_sd = std(Tf.db_cm); r.n = height(Tf);
            r.specimen = {sprintf('mean ± SD of %d', height(Tf))};
        else
            r.db_cm = Tf.db_cm; r.db_cm_sd = Tf.db_cm_sd; r.n = 1;
        end
        P = [P; r]; %#ok<AGROW>
    end
end
P.panel = repmat('b', height(P), 1);
P.panel(contains(P.layer, {'whole', 'temporal'})) = 'a';

%% style

S = paper_style();
studies = unique(P.study, 'stable');
cmap = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19; ...
        0.30 0.75 0.93; 0.64 0.08 0.18; 0.25 0.25 0.25; 0.85 0.45 0.70; 0.55 0.40 0.20; 0.00 0.60 0.50];
assert(numel(studies) <= size(cmap, 1), 'more studies than colours');
col = containers.Map(studies, num2cell(cmap(1:numel(studies), :), 2));
bandcol = [0.93 0.93 0.93];

%% plot

h = paper_figure([18 12]);
pos = {[1.3 5.4 6.6 6.0], [10.3 5.4 6.6 6.0]};
ttl = {'whole skull or layer not separated', 'cortical tables and diploë'};
for ip = 1:2
    ax = axes('Units', 'centimeters', 'Position', pos{ip}); hold(ax, 'on');
    pnl = char('a' + ip - 1);
    Pp = P(P.panel == pnl, :);
    top = Pp.db_cm + fillmissing(Pp.db_cm_sd, 'constant', 0);
    yl  = [0 1.05 * max([top(Pp.freq_mhz <= xl(2)); bnd(:)])];

    % background: neuromodulation band and simulation frequency
    patch(ax, fband([1 2 2 1]), [0 0 1e3 1e3], bandcol, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, [fsim fsim], [0 1e3], ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.5, 'HandleVisibility', 'off');

    hl = gobjects(0); lab = {};
    for s = unique(Pp.series, 'stable')'
        Q = sortrows(Pp(strcmp(Pp.series, s{1}), :), 'freq_mhz');
        c = col(Q.study{1});
        ls = '-';
        if contains(Q.layer{1}, 'diplo'); ls = '--'; end
        if Q.il(1); ls = ':'; end
        mfc = c; if Q.il(1); mfc = 'w'; end
        e = errorbar(ax, Q.freq_mhz, Q.db_cm, Q.db_cm_sd, 'LineStyle', ls, 'Color', c, ...
            'Marker', layer_marker(Q.layer{1}), 'MarkerSize', 3.5, 'MarkerFaceColor', mfc, ...
            'LineWidth', 0.8, 'CapSize', 2);
        hl(end+1) = e; %#ok<AGROW>
        lab{end+1} = Q.label{1}; %#ok<AGROW>
    end

    % Table 1 bounds at the simulation frequency (on top of the data)
    bw = 0.09;                                              % half width of the bars (MHz)
    hb = plot(ax, fsim + [-bw bw; -bw bw]', [bnd([1 3]); bnd([1 3])], '-', 'Color', 'k', 'LineWidth', 1.6);
    hr = plot(ax, fsim + [-bw bw], bnd([2 2]), '--', 'Color', 'k', 'LineWidth', 1.0);
    hl(end+1) = hb(1); lab{end+1} = 'Table 1 LOW / HIGH (0.5 MHz)'; %#ok<AGROW>
    hl(end+1) = hr;    lab{end+1} = 'Table 1 reference (0.5 MHz)';  %#ok<AGROW>

    xlim(ax, xl); ylim(ax, yl);
    set(ax, 'Layer', 'top');                                % axes above the shaded band
    xlabel(ax, 'frequency (MHz)'); ylabel(ax, 'attenuation (dB/cm)');
    title(ax, sprintf('(%s) %s', pnl, ttl{ip}));
    lg = legend(ax, hl, lab, 'Box', 'off', 'FontSize', 6, 'Interpreter', 'none', ...
        'Units', 'centimeters', 'NumColumns', 2);
    lg.ItemTokenSize = [14 6];
    lg.Position(1:2) = [pos{ip}(1) - 1.0, pos{ip}(2) - 1.3 - lg.Position(4)];   % below the axes
end

paper_save(h, pn.figures, 'suppfig3_skull_attenuation_frequency');
writetable(P(:, {'study', 'ref_key', 'label', 'specimen', 'layer', 'panel', 'il', 'freq_mhz', 'db_cm', 'db_cm_sd', 'n'}), ...
    fullfile(pn.figures, 'suppfig3_skull_attenuation_frequency_plotted.csv'));

%% local functions

function m = layer_marker(layer)
% marker per bone layer: outer table o, inner table s, cortical d, diploë ^
if contains(layer, 'outer');        m = 'o';
elseif contains(layer, 'inner');    m = 's';
elseif contains(layer, 'diplo');    m = '^';
elseif contains(layer, 'cortical'); m = 'd';
else;                               m = 'o';
end
end
