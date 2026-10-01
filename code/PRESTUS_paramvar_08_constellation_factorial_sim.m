clear all; close all; clc;

% PRESTUS_PARAMVAR_08_CONSTELLATION_FACTORIAL_SIM
%
% Benchmark 2 (brain+skull) ONLY. Computes the pairwise crossover between
% skull density and each of the other 13 tissue-property factors,
% motivated by PRESTUS_paramvar_09_morris_interactions_sim.m's pilot run
% flagging skull density as the only factor with a non-zero sigma
% (interaction signal) - based on just n=2 trajectories, too thin to trust
% on its own, hence this small, targeted, controlled follow-up.
%
% METHOD: for each of the other 13 factors (5 other acoustic-relevant:
% skull sound_speed/alpha_coeff, brain sound_speed/density/alpha_coeff;
% 8 thermal-relevant: thermal_conductivity/specific_heat_capacity/
% perfusion/absorption_fraction x brain/skull), a 2x2 grid is run: skull
% density (bottom/top) x that factor (bottom/top), with every OTHER
% property held at its literature midpoint. This is a direct, exact
% (not screened/statistical) interaction estimate per pair:
%   interaction = [outcome(density=top,factor=top)  - outcome(density=bottom,factor=top)]
%               - [outcome(density=top,factor=bottom) - outcome(density=bottom,factor=bottom)]
% If this is ~0, density's effect does not depend on that factor; if
% large, it does. 13 pairs x 4 combos = 52 evaluations total.
%
% NO ACOUSTIC/THERMAL CACHING (unlike PRESTUS_paramvar_09 and this
% script's own earlier version): every evaluation runs a full, independent
% fresh acoustic + thermal simulation. This was NOT the original design -
% an acoustic-cache-reuse scheme (2 acoustic "anchors" per density level,
% reused across the 8 thermal-factor pairs) was tried first and hit a real
% PRESTUS toolbox limitation: axisymmetric grids are expanded from their
% native 2D solve-space shape to a 3D analysis shape (convert_axisymmetric_
% to_3d) ONLY inside acoustic_wrapper, i.e. only on a fresh acoustic run.
% The cache-load path (cache_acoustic.m) never touches kwave_medium at
% all. Since these designs deliberately rebuild the medium fresh every
% evaluation (io.preproc_affix left unset, so new tissue-property values
% are actually picked up), a cache-reusing (acoustic-skipped) evaluation
% ends up combining a freshly-rebuilt, pre-expansion-shaped kwave_medium
% with an already-3D-expanded cached pressure field in acoustic_analysis -
% "Arrays have incompatible sizes for this operation". This is a toolbox-
% level gap (not something to patch there), specific to phantom +
% axisymmetric + heating + varying-tissue-properties-between-reuses, so it
% is worked around here by simply not attempting that reuse: every one of
% the 52 evaluations is a complete, independent simulation. This costs
% ~32 extra acoustic solves versus the cached design (still only 52 jobs
% total, well under budget), in exchange for correctness and simplicity -
% no phase ordering, no waiting, no acoustic-cache bookkeeping.
%
% JOB SUBMISSION: every evaluation is its own job, "auto mode" via
% prestus_pipeline_start(parameters) with parameters.platform set to
% hpc_detect_system()'s result ('slurm'/'qsub' if a scheduler is on PATH,
% else 'matlab' to run locally/synchronously) - same mechanism as every
% other script in this project. All 52 jobs are submitted directly; none
% depends on any other, so there is nothing to wait for.
%
% STORAGE: modules.run_nifti_creation / generate_report disabled;
% io.save_thermal_matrices / io.save_acoustic_matrices disabled throughout
% (nothing reuses any evaluation's cache, so nothing needs to persist it).
%
% SIMULATION DURATION: shortened thermal protocol (ptrd=20s, no post-
% stimulation cool-off) - see doc/PRESTUS_itrusst_methods.md, Analysis 4.
%
% ALWAYS validate with pilot_mode = true first (restricts to 2 of the 13
% factors = 8 evaluations) before committing the full 52-evaluation design.
%
% Output: data/tussim/<variant_name>/sub-002/sub-002_phantom<affix>.csv
% per evaluation.

user = 'julkos';

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
rootpath = fileparts(pathstr);

addpath(fullfile(rootpath, 'code'));
pn = setup_itrusst_paths(rootpath, user);

%% ── Configuration ───────────────────────────────────────────────────────────

subject_id   = 2;      % Benchmark 2 (brain+skull) only
setup_name   = 'itrusst_protocol1';
variant_name = 'itrusst_protocol1a_density_crossover';
pilot_mode   = false;   % true: restrict to a small validation subset (see header NOTE)

tissue_labels = {'brain', 'skull'};
tissue_suffixes = {'_b', '_s'};
acoustic_properties = {'sound_speed', 'density', 'alpha_coeff'};
acoustic_suffixes   = {'_c_', '_rho_', '_a0_'};
thermal_properties  = {'thermal_conductivity', 'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
thermal_suffixes    = {'_k_', '_heatc_', '_perf_', '_abs_'};

%% ── Detect job-submission platform (auto mode) ──────────────────────────────

try
    hpc_type = hpc_detect_system(); % 'slurm' or 'qsub' if a scheduler is on PATH
catch
    hpc_type = 'matlab';            % no scheduler found - run locally/synchronously
end
fprintf('Job submission platform: %s\n', hpc_type);

%% ── Load literature bounds + base config ────────────────────────────────────

params_table = readtable(fullfile(pn.data, 'tissue_parameters.csv'));

cd(fullfile(rootpath, 'data', 'configs'))
base_parameters = load_parameters(['config_', setup_name, '.yaml']);

if ismac
    base_parameters.startup.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
else
    base_parameters.startup.simnibs_bin_path = pn.simnibs;
end
pn.path_sim = fullfile(rootpath, 'data', 'tussim', variant_name);
if ~exist(pn.path_sim, 'dir'); mkdir(pn.path_sim); end
pn.figures = fullfile(rootpath, 'figures', 'paramvar');
if ~exist(pn.figures, 'dir'); mkdir(pn.figures); end
base_parameters.hpc.ld_library_path = "/opt/gcc/7.2.0/lib64";
base_parameters.path.anat = pn.path_seg;
base_parameters.path.seg  = pn.path_seg;
base_parameters.path.t1_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
base_parameters.path.t2_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
base_parameters.path.sim = pn.path_sim;
base_parameters.startup.paths_to_add = {pn.kwave, pn.minimize};
base_parameters.pct.enabled = 0;

base_parameters.platform = hpc_type;   % 'auto mode': whatever was detected above, resolved once here
base_parameters.simulation.medium = 'phantom';
base_parameters.simulation.interactive = 0;
base_parameters.simulation.code_type = 'cpp_cpu';
base_parameters.simulation.precision = 'single';
base_parameters.io.overwrite_files = 'never';
base_parameters.io.overwrite_simnibs = 0;
base_parameters.io.save_heatingvideo = 0;
base_parameters.io.save_matrices = 0;
base_parameters.io.save_thermal_matrices = 0;
base_parameters.io.save_acoustic_matrices = 0;
base_parameters.modules.run_nifti_creation = 0;
base_parameters.modules.generate_report = 0;
base_parameters.modules.run_source_setup  = 1;
base_parameters.modules.run_acoustic_sims = 1;
base_parameters.modules.run_heating_sims  = 1;
base_parameters.grid.resolution_mm = 0.5;
base_parameters.grid.default_dims = [70, 120] / base_parameters.grid.resolution_mm;
base_parameters.grid.pml_size = [20] / base_parameters.grid.resolution_mm;
base_parameters.grid.axisymmetric = 1;
base_parameters.subject_id = subject_id;

% Shortened thermal protocol (see header NOTE): 20s stimulation, 0s
% post-stimulation cool-off, instead of itrusst_protocol1's default 80s/80s.
base_parameters.timing.ptrd = 20;
base_parameters.timing.post_ptri_dur = 0;

base_parameters.hpc.timelimit = '00:30:00';  % see PRESTUS_paramvar_09's pilot run: 00:15:00 was NOT a
                                              % safe margin for the thermal stage at this protocol

%% ── Factor lists; skull density is the focal factor, all others are "other" ──

acoustic_factors = build_factor_list(tissue_labels, tissue_suffixes, acoustic_properties, acoustic_suffixes);
thermal_factors  = build_factor_list(tissue_labels, tissue_suffixes, thermal_properties, thermal_suffixes);

is_focal = strcmp({acoustic_factors.tissue}, 'skull') & strcmp({acoustic_factors.property}, 'density');
focal_factor = acoustic_factors(is_focal);
other_acoustic_factors = acoustic_factors(~is_focal);

levels = [1 2]; % 1 = bottom, 2 = top

%% ── Build the 52 full-evaluation job specs ─────────────────────────────────

specs = struct('eval_type', {}, 'other_factor', {}, 'extra_factors', {}, 'extra_levels', {}, 'output_affix', {});

% 20 acoustic-pair combos: density x each of the 5 other acoustic factors
for oa = 1:numel(other_acoustic_factors)
    other_factor = other_acoustic_factors(oa);
    for fl = levels
        for ol = levels
            pair_factors = [focal_factor, other_factor];
            pair_levels  = [fl, ol];
            specs(end+1) = struct('eval_type', "acoustic_pair", 'other_factor', other_factor, ...
                'extra_factors', pair_factors, 'extra_levels', pair_levels, ...
                'output_affix', combo_affix(pair_factors, pair_levels)); %#ok<AGROW>
        end
    end
end

% 32 thermal-pair combos: density x each of the 8 thermal factors
for tf = 1:numel(thermal_factors)
    thermal_factor = thermal_factors(tf);
    for fl = levels
        for ol = levels
            pair_factors = [focal_factor, thermal_factor];
            pair_levels  = [fl, ol];
            specs(end+1) = struct('eval_type', "thermal_pair", 'other_factor', thermal_factor, ...
                'extra_factors', pair_factors, 'extra_levels', pair_levels, ...
                'output_affix', combo_affix(pair_factors, pair_levels)); %#ok<AGROW>
        end
    end
end

if pilot_mode
    specs = specs([1 2 3 4 21 22 23 24]);  % first acoustic-pair factor's 4 combos + first thermal-pair factor's 4 combos
end

fprintf('Density crossover: %d evaluations%s.\n', numel(specs), ternary(pilot_mode, ' [PILOT MODE]', ''));

%% ── Design table (for PRESTUS_paramvar_13_density_crossover_plot.m) ─────────
% One row per evaluation. eval_type distinguishes acoustic_pair (density x
% one of the 5 other acoustic factors) from thermal_pair (density x one of
% the 8 thermal factors). other_tissue/other_property identify the
% non-focal factor in each pair.

n_rows = numel(specs);
design_rows = table('Size', [n_rows 6], ...
    'VariableTypes', {'string', 'string', 'string', 'double', 'double', 'string'}, ...
    'VariableNames', {'eval_type', 'other_tissue', 'other_property', 'focal_level', 'other_level', 'output_affix'});
for i = 1:n_rows
    other = specs(i).other_factor;
    design_rows(i, :) = {specs(i).eval_type, string(other.tissue), string(other.property), ...
        specs(i).extra_levels(1), specs(i).extra_levels(2), string(specs(i).output_affix)};
end

writetable(design_rows, fullfile(pn.figures, 'density_crossover_design.csv'));

%% ── Submit every evaluation (all independent - no phase ordering needed) ───

for i = 1:numel(specs)
    parameters = build_parameters(base_parameters, tissue_labels, acoustic_properties, ...
        thermal_properties, params_table, specs(i));
    prestus_pipeline_start(parameters);
end

%% ── Local functions ──────────────────────────────────────────────────────

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function factors = build_factor_list(tissue_labels, tissue_suffixes, property_labels, property_suffixes)
    factors = struct('tissue', {}, 'tissue_suffix', {}, 'property', {}, 'property_suffix', {});
    for t = 1:numel(tissue_labels)
        for p = 1:numel(property_labels)
            factors(end+1) = struct('tissue', tissue_labels{t}, 'tissue_suffix', tissue_suffixes{t}, ...
                'property', property_labels{p}, 'property_suffix', property_suffixes{p}); %#ok<AGROW>
        end
    end
end

function affix = combo_affix(factors, levels)
    affix = '';
    for f = 1:numel(factors)
        affix = [affix factors(f).tissue_suffix factors(f).property_suffix num2str(levels(f))]; %#ok<AGROW>
    end
end

function parameters = set_mid_properties(parameters, tissue_labels, property_labels, params_table)
    for t = 1:numel(tissue_labels)
        for p = 1:numel(property_labels)
            parameters.medium_properties.(tissue_labels{t}).(property_labels{p}) = ...
                params_table.(property_labels{p})(strcmp(params_table.tissue, tissue_labels{t}) & ...
                                                    strcmp(params_table.level, 'mid'));
        end
    end
end

function parameters = apply_factor_levels(parameters, factors, levels, params_table)
    level_names = {'bottom', 'top'};
    for f = 1:numel(factors)
        parameters.medium_properties.(factors(f).tissue).(factors(f).property) = ...
            params_table.(factors(f).property)(strcmp(params_table.tissue, factors(f).tissue) & ...
                                                 strcmp(params_table.level, level_names{levels(f)}));
    end
end

function parameters = build_parameters(base_parameters, tissue_labels, acoustic_properties, ...
    thermal_properties, params_table, spec)
% Full, independent evaluation: density x one other factor, everything
% else at midpoint. Always a fresh acoustic + thermal solve (see header
% NOTE on why acoustic caching is not used here).
    parameters = base_parameters;
    parameters = set_mid_properties(parameters, tissue_labels, [acoustic_properties, thermal_properties], params_table);
    parameters = apply_factor_levels(parameters, spec.extra_factors, spec.extra_levels, params_table);
    parameters.io.output_affix = spec.output_affix;
end
