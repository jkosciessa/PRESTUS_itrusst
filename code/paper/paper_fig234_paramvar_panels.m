% PAPER_FIG234_PARAMVAR_PANELS
%
% Range panels of the one-at-a-time parameter variation (itrusst_protocol1a)
% for paper Figs. 2-4, with dashed reference-medium lines, plus the
% BabelBrain counterparts of the sub-002 panels (benchmark 2, Sam Pichardo,
% data/pichardo; load_babelbrain_benchmark2) on a shared y-range.
%
%   Fig. 2  brain I_SPPA, acoustic properties of brain and skull (sub-002)
%   Fig. 3a skull temperature rise, skull properties (sub-002)
%   Fig. 3b temperature rise at the focus, skull properties (sub-002)
%   Fig. 3c temperature rise at the focus, brain properties, without skull
%           (sub-001; red line: reference with skull) and with skull (sub-002)
%   Fig. 4  brain I_SPPA and skull temperature rise for the 2.5 mm skull
%           (sub-352, maximal transmission; red line: 2.5 mm reference,
%           grey: 6.5 mm reference)
%
% Reference lines come from the reference-medium runs (paper_reference_metrics);
% thermal references need a baseline run (PRESTUS_paramvar_00_baseline_sim.m),
% which currently exists for sub-002 only.
%
% Focal heating is read at the empirical acoustic focus (riseT_focus, see
% load_focus_heating_metrics); set focus_metric = 'riseT_target' for the
% fixed nominal target. BabelBrain reports the temperature rise at the
% target only ('Max. Temperature Target'), which its focal panels show
% regardless of focus_metric.
%
% BabelBrain was run at 6 and 9 points per wavelength (PPW); bb_ppw selects
% the run shown. 9 PPW is BabelBrain's best available estimate; 6 PPW
% (~ the PRESTUS 0.5 mm grid in water at 500 kHz) serves as convergence
% check: intensities change by <= 3 %, peak heating by up to ~30 %.
% Differences between the two runs and PRESTUS:
% PRESTUS_paramvar_02_benchmark2_babelbrain_plot.m.
%
% Output: figures/paper/fig{2,3,4}/*.{pdf,png}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn.tussim   = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a');
pn.pichardo = itrusst_input_path(fullfile(rootpath, 'data', 'pichardo'), ...
    'ITRUSST_BABELBRAIN_DIR', 'BabelBrain benchmark 2 results (data/pichardo)');
pn.figures  = fullfile(rootpath, 'figures', 'paper');

%% settings

focus_metric = 'riseT_focus';
bb_ppw       = 9;

prop.labels   = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                 'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
prop.suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};
prop.symbols  = {'c', '\rho', '\alpha', '\kappa', 'c_p', '\omega', '\beta'};
acoustic      = 1:3;
sweep         = 1:2; % bottom, top

S = paper_style();

%% load

dir_sub = @(id) fullfile(pn.tussim, sprintf('sub-%03d', id));
load_res = @(id, tissue, idx) load_paramvar_results(dir_sub(id), id, {tissue}, {['_' tissue(1)]}, ...
    prop.labels(idx), prop.suffixes(idx), sweep, prop.symbols(idx));
load_focus = @(id, tissue, idx) load_focus_heating_metrics(dir_sub(id), id, {tissue}, {['_' tissue(1)]}, ...
    prop.labels(idx), prop.suffixes(idx), sweep, prop.symbols(idx));

ref001 = paper_reference_metrics(dir_sub(1), 1);
ref002 = paper_reference_metrics(dir_sub(2), 2);
ref352 = paper_reference_metrics(dir_sub(352), 352);

b002a = load_res(2, 'brain', acoustic);  s002a = load_res(2, 'skull', acoustic);
s002  = load_res(2, 'skull', 1:7);
f002s = load_focus(2, 'skull', 1:7);     f002b = load_focus(2, 'brain', 1:7);
f001b = load_focus(1, 'brain', 1:7);
b352a = load_res(352, 'brain', acoustic); s352a = load_res(352, 'skull', acoustic);
s352  = load_res(352, 'skull', 1:7);

bb_load = @(tissue, idx) load_babelbrain_benchmark2(pn.pichardo, bb_ppw, {tissue}, prop.labels(idx));
[bb_b_a, bb_ref] = bb_load('brain', acoustic);  bb_s_a = bb_load('skull', acoustic);
bb_s = bb_load('skull', 1:7);                   bb_b   = bb_load('brain', 1:7);

ref = @(v, c, varargin) struct('value', v, 'color', c, 'label', {strjoin(varargin, '')});
bb  = @(data, refs) struct('data', data, 'refs', refs);
sym_ac = prop.symbols(acoustic);
tis_ac = [repmat({'b'}, 1, 3), repmat({'s'}, 1, 3)];

%% panels: {figure, name, data, symbols, tissues, ylabel, ylim, refs, BabelBrain}
% BabelBrain: struct(data, refs) of the counterpart panel, or [] if none

panels = {
  'fig2', 'fig2_prestus_isppa_brain', [b002a.Isppa_brain, s002a.Isppa_brain], [sym_ac sym_ac], tis_ac, ...
          'Brain I_{SPPA} (W/cm^2)', [], ref(ref002.Isppa_brain, S.ref), ...
          bb([bb_b_a.Isppa_brain, bb_s_a.Isppa_brain], ref(bb_ref.Isppa_brain, S.ref))
  'fig3', 'fig3a_prestus_riseT_skull', s002.riseT_skull, prop.symbols, repmat({'s'}, 1, 7), ...
          'Skull temperature rise (°C)', [], ref(ref002.riseT_skull, S.ref), ...
          bb(bb_s.riseT_skull, ref(bb_ref.riseT_skull, S.ref))
  'fig3', 'fig3b_prestus_riseT_focus_skullprops', f002s.(focus_metric), prop.symbols, repmat({'s'}, 1, 7), ...
          'Focal temperature rise (°C)', [], ref(ref002.(focus_metric), S.ref), ...
          bb(bb_s.riseT_target, ref(bb_ref.riseT_target, S.ref))
  'fig3', 'fig3c_prestus_riseT_focus_noskull', f001b.(focus_metric), prop.symbols, repmat({'b'}, 1, 7), ...
          'Focal temperature rise (°C)', [], [ref(ref001.(focus_metric), S.ref), ...
          ref(ref002.(focus_metric), S.withskull, 'with skull (reference)')], []
  'fig3', 'fig3c_prestus_riseT_focus_brainprops', f002b.(focus_metric), prop.symbols, repmat({'b'}, 1, 7), ...
          'Focal temperature rise (°C)', [], ref(ref002.(focus_metric), S.ref), ...
          bb(bb_b.riseT_target, ref(bb_ref.riseT_target, S.ref))
  'fig4', 'fig4b_prestus_isppa_brain_2p5mm', [b352a.Isppa_brain, s352a.Isppa_brain], [sym_ac sym_ac], tis_ac, ...
          'Brain I_{SPPA} (W/cm^2)', [], [ref(ref352.Isppa_brain, S.highlight), ...
          ref(ref002.Isppa_brain, S.ref_alt, '6.5 mm default')], []
  'fig4', 'fig4c_prestus_riseT_skull_2p5mm', s352.riseT_skull, prop.symbols, repmat({'s'}, 1, 7), ...
          'Skull temperature rise (°C)', [], [ref(ref352.riseT_skull, S.highlight), ...
          ref(ref002.riseT_skull, S.ref_alt, '6.5 mm default')], []
};

%% plot

for i = 1:size(panels, 1)
    [fig, name, data, syms, tis, ylab, ylims, refs, bbp] = panels{i, :};
    if all(isnan(data(:)))
        warning('%s: no data, skipped.', name);
        continue;
    end
    h = paper_figure(S.range_fig);
    ax = axes(h, 'Units', 'centimeters', 'Position', S.range_axes);
    paper_range_panel(ax, data, syms, tis, ylab, ylims, refs);
    if ~isempty(bbp)
        % BabelBrain counterpart; both panels share the larger y-range
        hb = paper_figure(S.range_fig);
        axb = axes(hb, 'Units', 'centimeters', 'Position', S.range_axes);
        paper_range_panel(axb, bbp.data, syms, tis, ylab, ylims, bbp.refs);
        if isempty(ylims)
            set([ax axb], 'YLim', [0, max([ax.YLim(2), axb.YLim(2)])]);
        end
        paper_save(hb, fullfile(pn.figures, fig), strrep(name, 'prestus', 'babelbrain'));
        close(hb);
    end
    paper_save(h, fullfile(pn.figures, fig), name);
    close(h);
end
