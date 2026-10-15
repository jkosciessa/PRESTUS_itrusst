function S = paper_style()
% Shared look of the paper panels (code/paper/): font, colours, sizes.
% Panels are exported at print size (cm), so font sizes are final sizes.

S.font      = 'Helvetica';
S.fontsize  = 8;
S.linewidth = 1;

% tissues (cf. plot_metric_with_patches)
S.brain     = [0.50 0.50 0.50];
S.skull     = [0.70 0.00 0.00];
S.water     = [0.10 0.75 0.85];

% benchmarks: no skull (sub-001) / with skull (sub-002)
S.noskull   = [0.00 0.00 0.00];
S.withskull = [0.85 0.10 0.10];

% sweep levels: bottom (low) / top (high)
S.low       = [1 1 1];
S.high      = [0 0 0];

% low / high level as line colours (supplementary figures)
S.low_line  = [0.00 0.45 0.74];
S.high_line = [0.85 0.10 0.10];

% reference lines
S.ref       = [0.00 0.00 0.00];   % reference medium, same benchmark
S.ref_alt   = [0.55 0.55 0.55];   % reference of another benchmark
S.highlight = [0.90 0.10 0.10];   % 2.5 mm (maximal transmission) benchmark skull

% setups and toolbox defaults (Fig. 5)
S.setup     = [0.10 0.10 0.10];
S.PRESTUS   = [0.15 0.50 0.80];
S.kPlan     = [0.75 0.22 0.22];
S.BabelBrain= [0.05 0.55 0.25];
S.bar_int   = [0.60 0.60 0.60];
S.bar_temp  = [0.70 0.20 0.20];

% panel geometry [left bottom width height] of the axes box in cm, so that
% panels of the same kind line up when placed side by side
S.range_axes   = [1.3 1.2 3.4 4.2];  % paramvar range panels
S.range_fig    = [5.0 5.8];          % figure size (cm)
S.profile_axes = [1.0 1.0 2.4 4.6];  % vertical axial profiles
S.profile_fig  = [4.2 6.3];
end
