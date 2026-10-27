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
% The yellow band marks the skull attenuation LOW-HIGH range of Table 1 at
% 0.5 MHz (dashed: reference), read from data/tissue_parameters.csv and drawn
% across the width of the grey band. Both panels share the y axis. In (a),
% series of one study share a hue (light to dark shades); in (b), colour marks
% the layer (outer table, cortical with tables not separated, inner table:
% dark to light blue; diploë: orange) and the marker the study. Legend order
% is chronological.
%
% Two versions: 0-2.3 MHz (all values; lines run on beyond the axis), and
% 0-1 MHz with the values above 1 MHz left out (series without values up to
% 1 MHz drop out; the script lists them).
%
% Output: figures/paper/supp_fig3/suppfig3_skull_attenuation_frequency{,_0-1MHz}.{pdf,png}
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

% series in chronological order of the studies (legend order follows)
yr = cellfun(@(s) str2double(regexp(s, '\d{4}', 'match', 'once')), P.study);
[~, io] = sort(yr);                                         % sort is stable
P = P(io, :);

% (a) one hue per study (chronological); several series of one study get
% shades of that hue, from lighter to darker
studies = unique(P.study, 'stable');
cmap = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.00 0.60 0.50; 0.49 0.18 0.56; 0.64 0.08 0.18; ...
        0.25 0.25 0.25; 0.55 0.40 0.20; 0.85 0.33 0.10; 0.30 0.75 0.93; 0.85 0.45 0.70; 0.93 0.69 0.13];  % no yellow before the 11th study: yellow marks Table 1
assert(numel(studies) <= size(cmap, 1), 'more studies than colours');
col = containers.Map(studies, num2cell(cmap(1:numel(studies), :), 2));
shade = [-0.35 0.5];                                        % lightest (mix with white) / darkest (mix with black)

% (b) colour by layer (cortical blues, diploë orange), marker by study;
% several series of one study and layer get shades of the layer colour
lcol = struct('outer', [0.03 0.27 0.58], 'cortical', [0.13 0.47 0.80], ...
              'inner', [0.42 0.68 0.90], 'diploe', [0.90 0.40 0.05]);
shade_b = [-0.3 0.35];
smark = {'o', 'd', 's', '^', 'v', 'p', 'h'};                % markers per study in (b), chronological
bstud = unique(P.study(P.panel == 'b'), 'stable');
assert(numel(bstud) <= numel(smark), 'more studies in (b) than markers');
mrk = containers.Map(bstud, smark(1:numel(bstud)));

bandcol = [0.93 0.93 0.93];
bndcol  = [1.00 0.92 0.45];                                 % Table 1 LOW-HIGH range

%% plot: full range, and 0-1 MHz with the values above 1 MHz left out

V = struct('fmax', {2.3, 1}, 'xtick', {0:0.5:2, 0:0.2:1}, 'crop', {false, true}, ...
           'name', {'suppfig3_skull_attenuation_frequency', 'suppfig3_skull_attenuation_frequency_0-1MHz'});
for iv = 1:numel(V)
    xl = [0 V(iv).fmax];
    Pv = P;
    if V(iv).crop; Pv = P(P.freq_mhz <= V(iv).fmax, :); end

    h = paper_figure([18 12]);
    x0  = [1.3 10.3];                                       % left edge of the axes (cm)
    w   = 6.9; ah = 5.6;                                    % axes width / height (cm)
    ttl = {'whole skull or layer not separated', 'cortical tables (blue) and diploë (orange)'};
    ax = gobjects(1, 2); lg = gobjects(1, 2);
    ymax = max(bnd);
    for ip = 1:2
        ax(ip) = axes('Units', 'centimeters', 'Position', [x0(ip) 1 w ah]); hold(ax(ip), 'on');
        pnl = char('a' + ip - 1);
        Pp = Pv(Pv.panel == pnl, :);
        top = Pp.db_cm + fillmissing(Pp.db_cm_sd, 'constant', 0);
        ymax = max([ymax; top(Pp.freq_mhz <= xl(2))]);

        % background: neuromodulation band, Table 1 skull attenuation at the
        % simulation frequency (LOW-HIGH shaded across the band, reference
        % dashed) and the simulation frequency, all beneath the data
        patch(ax(ip), fband([1 2 2 1]), [0 0 1e3 1e3], bandcol, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        patch(ax(ip), fband([1 2 2 1]), bnd([1 1 3 3]), bndcol, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax(ip), fband, bnd([2 2]), '--', 'Color', 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        plot(ax(ip), [fsim fsim], [0 1e3], ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.5, 'HandleVisibility', 'off');

        series = unique(Pp.series, 'stable');
        first  = cellfun(@(s) find(strcmp(Pp.series, s), 1), series);
        grp    = Pp.study(first);                           % colour group: study (a), study and layer (b)
        if pnl == 'b'; grp = strcat(grp, '|', cellfun(@layer_class, Pp.layer(first), 'UniformOutput', false)); end
        hl = gobjects(0); lab = {};
        for is = 1:numel(series)
            Q = sortrows(Pp(strcmp(Pp.series, series{is}), :), 'freq_mhz');
            if pnl == 'a'
                c = col(Q.study{1}); sh = shade; m = layer_marker(Q.layer{1});
            else
                c = lcol.(layer_class(Q.layer{1})); sh = shade_b; m = mrk(Q.study{1});
            end
            same = find(strcmp(grp, grp{is}));
            if numel(same) > 1
                t = sh(1) + diff(sh) * (find(same == is) - 1) / (numel(same) - 1);
                if t < 0; c = c + (1 - c) * -t; else; c = c * (1 - t); end
            end
            ls = '-';
            if contains(Q.layer{1}, 'diplo'); ls = '--'; end
            if Q.il(1); ls = ':'; end
            mfc = c; if Q.il(1); mfc = 'w'; end
            hl(end+1) = errorbar(ax(ip), Q.freq_mhz, Q.db_cm, Q.db_cm_sd, 'LineStyle', ls, 'Color', c, ...
                'Marker', m, 'MarkerSize', 3.5, 'MarkerFaceColor', mfc, 'LineWidth', 0.8, 'CapSize', 2); %#ok<AGROW>
            lab{end+1} = Q.label{1}; %#ok<AGROW>
        end

        set(ax(ip), 'Layer', 'top', 'XTick', V(iv).xtick);  % axes above the shaded band
        xlim(ax(ip), xl);
        xlabel(ax(ip), 'frequency (MHz)'); ylabel(ax(ip), 'attenuation (dB/cm)');
        title(ax(ip), sprintf('(%s) %s', pnl, ttl{ip}));
        lg(ip) = legend(ax(ip), hl, lab, 'Box', 'off', 'FontSize', 6, 'Interpreter', 'none', ...
            'Units', 'centimeters', 'NumColumns', 1);
        lg(ip).ItemTokenSize = [14 6];
    end
    set(ax, 'YLim', [0 1.05 * ymax]);                       % common y axis, panels comparable

    % key shared by both panels (dummy handles in a hidden axes)
    axk = axes('Units', 'centimeters', 'Position', [0 0 1 1], 'Visible', 'off'); hold(axk, 'on');
    gk = [0.45 0.45 0.45];
    hk = [patch(axk, NaN, NaN, bandcol, 'EdgeColor', 'none'), ...
          plot(axk, NaN, NaN, '--', 'Color', gk, 'LineWidth', 0.8), ...
          patch(axk, NaN, NaN, bndcol, 'EdgeColor', 'none'), ...
          plot(axk, NaN, NaN, '--', 'Color', 'k', 'LineWidth', 0.8), ...
          plot(axk, NaN, NaN, ':o',  'Color', gk, 'MarkerFaceColor', 'w', 'MarkerSize', 3.5, 'LineWidth', 0.8)];
    labk = {sprintf('common TUS frequencies (%.1f–%.2f MHz)', fband), ...
            'dashed: diploë', ...
            sprintf('Table 1 skull LOW–HIGH (at %.1f MHz)', fsim), ...
            sprintf('Table 1 skull reference (at %.1f MHz)', fsim), ...
            'open, dotted: includes reflection [IL]'};
    lgk = legend(axk, hk, labk, 'Box', 'off', 'FontSize', 6, 'Units', 'centimeters', 'NumColumns', 3);
    lgk.ItemTokenSize = [22 6];                             % long enough to show the dashes

    % layout from the bottom: series legends, axes, titles, key
    lh = max(arrayfun(@(l) l.Position(4), lg));
    yb = 0.2 + lh + 1.2;                                    % axes bottom, room for the x label
    for ip = 1:2
        ax(ip).Position = [x0(ip) yb w ah];
        lg(ip).Position(1:2) = [x0(ip) - 1.0, yb - 1.1 - lg(ip).Position(4)];
    end
    yk = yb + ah + 0.7;                                     % key above the titles
    lgk.Position(1:2) = [(18 - lgk.Position(3)) / 2, yk];
    h.Position(4) = yk + lgk.Position(4) + 0.1;

    paper_save(h, pn.figures, V(iv).name);
    writetable(Pv(:, {'study', 'ref_key', 'label', 'specimen', 'layer', 'panel', 'il', 'freq_mhz', 'db_cm', 'db_cm_sd', 'n'}), ...
        fullfile(pn.figures, [V(iv).name '_plotted.csv']));
    dropped = setdiff(unique(P.label), unique(Pv.label));
    if ~isempty(dropped); fprintf('  %s: no values up to %g MHz for %s\n', V(iv).name, V(iv).fmax, strjoin(dropped, '; ')); end
end

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

function c = layer_class(layer)
% colour class per bone layer in (b): outer, inner, cortical (tables not separated), diploe
if contains(layer, 'outer');        c = 'outer';
elseif contains(layer, 'inner');    c = 'inner';
elseif contains(layer, 'diplo');    c = 'diploe';
else;                               c = 'cortical';
end
end
