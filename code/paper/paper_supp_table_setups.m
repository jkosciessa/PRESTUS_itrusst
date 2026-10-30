% PAPER_SUPP_TABLE_SETUPS
%
% Supplementary table: brain and skull properties of the literature setups
% (Temp-H / Temp-L = heat_max / heat_min, Trans-H / Trans-L = trans_max /
% trans_min) and the emulated toolbox defaults (PRESTUS, k-Plan*,
% BabelBrain*), as simulated in the pipeline comparison
% (PRESTUS_pipeline_01_benchmark1_sim.m; Fig. 5). The reference run is not
% shown in Fig. 5 and its values are in Table 1, so it is omitted. Trans
% setups bracket intensity only: absorption and thermal rows are shown as '–'
% (held at the reference).
%
% Values are read from data/configs/config_medium_<setup>.yaml, i.e. the media
% the simulations load (heat_*/trans_* are generated from
% data/tissue_parameters.csv + tissue_setups.csv, with brain absorption as a
% fraction of attenuation and skull absorption capped at attenuation; the tool
% configs are hand-written). Attenuation and absorption
% are given at 500 kHz in dB/cm (alpha_coeff * 0.5^alpha_power); beta/alpha
% is the absorbed fraction of attenuation implied by the two coefficients.
% k-Plan perfusion is not reported and simulated as 0 (shown as n.r.).
% Water is identical across setups and given in a footnote.
%
% Output: figures/paper/supp_table/supptable_setups.{csv,md}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath); % PRESTUS yaml reader
pn.out = fullfile(rootpath, 'figures', 'paper', 'supp_table');
if ~exist(pn.out, 'dir'); mkdir(pn.out); end

%% setup properties

freq_mhz = 0.5;
setups = {'heat_max', 'heat_min', 'trans_max', 'trans_min', 'PRESTUS', 'kPlan', 'BabelBrain'};
labels = {'Temp-H', 'Temp-L', 'Trans-H', 'Trans-L', 'PRESTUS', 'k-Plan*', 'BabelBrain*'};
% rows not varied in the transmission setups (held at the reference)
heat_only = {'absorption_500', 'absorption_ratio', 'thermal_conductivity', ...
             'specific_heat_capacity', 'perfusion'};

media = cell(1, numel(setups));
for s = 1:numel(setups)
    cfg = yaml.loadFile(fullfile(pn.configs, ['config_medium_' setups{s} '.yaml']));
    media{s} = cfg.medium_properties;
end

% row: field (or derived quantity), label, unit
rows = {'sound_speed',            'Sound speed c',                      'm/s'; ...
        'density',                'Density ρ',                          'kg/m³'; ...
        'alpha_500',              'Attenuation α at 500 kHz',           'dB/cm'; ...
        'absorption_500',         'Absorption β at 500 kHz',            'dB/cm'; ...
        'absorption_ratio',       'Absorbed fraction β/α',              '–'; ...
        'alpha_power',            'Frequency power law y',              '–'; ...
        'thermal_conductivity',   'Thermal conductivity κ',             'W/m/K'; ...
        'specific_heat_capacity', 'Specific heat capacity c_p',         'J/kg/K'; ...
        'perfusion',              'Perfusion ω',                        'mL/min/kg'};

tissues = {'brain', 'Brain'; 'skull', 'Skull'};

cells = {};
for t = 1:size(tissues, 1)
    for r = 1:size(rows, 1)
        line = [tissues(t, 2), rows(r, 2), rows(r, 3)];
        for s = 1:numel(setups)
            m = media{s}.(tissues{t, 1});
            switch rows{r, 1}
                case 'alpha_500'
                    v = m.alpha_coeff * freq_mhz ^ m.alpha_power;
                case 'absorption_500'
                    v = m.absorption_coeff * freq_mhz ^ m.alpha_power;
                case 'absorption_ratio'
                    v = m.absorption_coeff / m.alpha_coeff;
                otherwise
                    v = m.(rows{r, 1});
            end
            if startsWith(setups{s}, 'trans') && ismember(rows{r, 1}, heat_only)
                line{end+1} = '–'; %#ok<SAGROW> held at the reference
            elseif strcmp(setups{s}, 'kPlan') && strcmp(rows{r, 1}, 'perfusion')
                line{end+1} = 'n.r.'; %#ok<SAGROW> not reported, simulated as 0
            elseif strcmp(rows{r, 1}, 'absorption_ratio')
                line{end+1} = sprintf('%.2f', v); %#ok<SAGROW>
            else
                line{end+1} = format_value(v); %#ok<SAGROW>
            end
        end
        cells(end+1, :) = line; %#ok<SAGROW>
    end
end
header = [{'Tissue', 'Property', 'Unit'}, labels];

%% footnotes

w = media{1}.water;
for s = 2:numel(setups)
    if ~isequal(media{s}.water, w)
        error('Water differs between %s and %s.', setups{1}, setups{s});
    end
end
notes = { ...
    ['Attenuation and absorption at 500 kHz; PRESTUS stores the 1 MHz coefficients ' ...
     'α₀ = α(500 kHz) / 0.5^y and β₀ = β(500 kHz) / 0.5^y in dB/(MHz^y·cm); attenuation and absorption share y. ' ...
     'Heat deposition uses β in all setups (thermal.absorption_model: coefficient).'], ...
    ['Temp-H / Temp-L: per property, the literature bound that maximises / minimises ' ...
     'predicted brain heating (data/tissue_setups.csv). Trans-H / Trans-L: per acoustic ' ...
     'property, the bound that maximises / minimises transcranial transmission; absorption ' ...
     'and thermal properties are held at the reference values (Table 1; –), as these setups ' ...
     'bracket intensity only. Brain absorption is a fraction ' ...
     'of the attenuation it is combined with (29 / 77 / 100 %); skull absorption is a ' ...
     'coefficient, capped at the attenuation.'], ...
    'PRESTUS: v05 defaults (config/config_default.yaml).', ...
    ['* Emulated in PRESTUS from the tool defaults. Both tools specify absorption as a ' ...
     'fraction of attenuation; β is that fraction times the tool''s own attenuation. ' ...
     'k-Plan skull sound speed follows its density mapping (1.33 ρ + 166.7) at the ' ...
     'thermal-simulation density of 1850 kg/m³. n.r.: not reported by k-Plan, simulated as 0.'], ...
    sprintf(['Water (all setups): c = %s m/s, ρ = %s kg/m³, α = %s dB/cm at 500 kHz ' ...
     '(y = %s), β = %s dB/cm, κ = %s W/m/K, c_p = %s J/kg/K, no perfusion.'], ...
     format_value(w.sound_speed), format_value(w.density), ...
     format_value(w.alpha_coeff * freq_mhz ^ w.alpha_power), format_value(w.alpha_power), ...
     format_value(w.absorption_coeff * freq_mhz ^ w.alpha_power), ...
     format_value(w.thermal_conductivity), format_value(w.specific_heat_capacity))};

%% write

writecell([header; cells], fullfile(pn.out, 'supptable_setups.csv'));

fid = fopen(fullfile(pn.out, 'supptable_setups.md'), 'w', 'n', 'UTF-8');
fprintf(fid, '| %s |\n', strjoin(header, ' | '));
fprintf(fid, '|%s\n', repmat('---|', 1, numel(header)));
for r = 1:size(cells, 1)
    fprintf(fid, '| %s |\n', strjoin(cells(r, :), ' | '));
end
fprintf(fid, '\n');
for n = 1:numel(notes)
    fprintf(fid, '%s\n\n', notes{n});
end
fclose(fid);

fprintf('Wrote %s\n', fullfile(pn.out, 'supptable_setups.{csv,md}'));

function str = format_value(v)
% integers above 100 without decimals, otherwise three significant digits
if abs(v) >= 100
    str = sprintf('%.0f', v);
else
    str = sprintf('%.3g', v);
end
end
