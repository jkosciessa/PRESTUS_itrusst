clear all; close all; clc;

% PRESTUS_PARAMVAR_09_MORRIS_INTERACTIONS_SIM
%
% EXPLORATORY (side analysis, not part of the main chain or
% itrusst_rerun_all). Benchmark 2 (brain+skull) ONLY. A Morris elementary-
% effects (one-factor-at-a-time trajectory) screening design across all 14
% tissue properties, to measure INTERACTION strength per property.
% PRESTUS_paramvar_08_constellation_factorial_sim.m gives exact pairwise
% interactions, but only for skull density x each other factor; this
% design screens every factor's context-dependence at once.
%
% METHOD (Morris 1991 elementary effects, on the same 2-level bottom/top
% literature-bound grid used throughout this project): r independent
% TRAJECTORIES are drawn. Each starts at a random point in the 14-factor
% {bottom,top} hypercube, then flips one not-yet-flipped factor at a time
% (each factor exactly once per trajectory, in a random order) until all
% 14 have been flipped - a "staircase" walk of 15 points (1 baseline + 14
% steps). At step s, the ELEMENTARY EFFECT of the flipped factor is
% outcome(point_s) - outcome(point_{s-1}), signed so it always reads as
% "effect of bottom -> top" regardless of which direction this trajectory
% happened to flip it. After r trajectories, each factor has r elementary
% effects:
%   mu*  = mean(|elementary effects|)  - overall importance (like Analysis
%          2's OAT effect, but now measured with every OTHER factor at a
%          different random point on every trajectory, not held at mid)
%   sigma = std(elementary effects)    - INTERACTION/context-dependence:
%          if a factor's effect is the same regardless of where the other
%          13 factors are, sigma ~ 0; if its effect depends on the other
%          factors' current values, sigma is large. This is the direct,
%          per-factor interaction measure this design was built to get.
% (With only 2 levels per factor there is no separate "curvature" term to
% confound sigma with, unlike continuous Morris designs with >=3 levels -
% sigma here is a clean interaction signal.)
%
% FACTOR ORDER: within each trajectory, the 6 acoustic-relevant factors
% (sound_speed/density/alpha_coeff x brain/skull) are flipped FIRST (in
% random relative order among themselves), then the 8 thermal-relevant
% factors (thermal_conductivity/specific_heat_capacity/perfusion/
% absorption_fraction x brain/skull) - a grouped ordering (a valid Morris
% design choice). It was introduced for an acoustic-cache-reuse scheme and
% is kept so that the design (rng_seed 42) stays the same as before.
%
% NO ACOUSTIC/THERMAL CACHING: every evaluation is a full, independent,
% fresh acoustic + thermal simulation (as in PRESTUS_paramvar_08). The
% earlier version of this script reused the acoustic cache of point 7 for
% the 8 thermal-only points of each trajectory; on the axisymmetric
% phantom grid every one of these cache-reusing jobs failed in
% acoustic_analysis ("Arrays have incompatible sizes", cached 3D-expanded
% pressure field vs. freshly built 2D medium; see PRESTUS_paramvar_08's
% header). That run (data/tussim_261001/itrusst_protocol1a_morris, old
% tissue ranges) therefore ended with 47/150 results and no elementary
% effects for any thermal factor.
%
% JOB SUBMISSION: every evaluation is its own job, submitted via
% prestus_pipeline_start(parameters) with parameters.platform set to
% whatever hpc_detect_system() finds ('slurm'/'qsub' if a scheduler is on
% PATH, else running locally/synchronously). No job depends on another, so
% all are submitted at once.
%
% COST: 15 evaluations per trajectory; r=10 (default) gives 150 jobs.
% Reference timings (Benchmark 2, 0.5 mm axisymmetric, shortened thermal
% protocol): 6-29 min per fresh evaluation (mean ~19 min in the
% PRESTUS_paramvar_08 run), i.e. ~30-50 CPU job-hours (cpp_cpu) for r=10. Increase r
% for more reliable mu*/sigma estimates (cost scales linearly).
%
% Points shared between trajectories (same 14-factor combination, hence the
% same output_affix) are submitted once.
%
% STORAGE: modules.run_nifti_creation / generate_report disabled;
% io.save_thermal_matrices / io.save_acoustic_matrices disabled (nothing
% reuses any evaluation's cache).
%
% ANALYSIS: PRESTUS_paramvar_10_morris_analysis_plot.m joins the design
% table this script writes (data/tussim/<variant_name>/morris_design.csv;
% columns: trajectory, step, tissue, property, direction, output_affix)
% against each output_affix's sub-002_phantom<affix>.csv.
%
% ALWAYS validate with pilot_mode = true (restricts to 2 trajectories, 30
% evaluations) before committing the full design.
%
% Output: data/tussim/<variant_name>/sub-002/sub-002_phantom<affix>.csv
% per evaluation; data/tussim/<variant_name>/morris_design.csv (trajectory design
% table for the follow-up analysis).

user = 'julkos';

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
rootpath = fileparts(fileparts(pathstr));

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, user);

%% ── Configuration ───────────────────────────────────────────────────────────

subject_id   = 2;      % Benchmark 2 (brain+skull) only
setup_name   = 'itrusst_protocol';
variant_name = 'itrusst_protocol1a_morris';
pilot_mode   = false;   % true: restrict to 2 trajectories (30 evaluations) for validation
r            = 10;     % number of Morris trajectories; cost scales as 15*r evaluations
rng_seed     = 42;     % fixed seed so the design is reproducible across reruns

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
base_parameters.modules.run_source_setup  = 1;
base_parameters.modules.run_acoustic_sims = 1;
base_parameters.modules.run_heating_sims  = 1;
base_parameters.modules.run_nifti_creation = 0;
base_parameters.modules.generate_report = 0;
base_parameters.grid.resolution_mm = 0.5;
base_parameters.grid.default_dims = [70, 120] / base_parameters.grid.resolution_mm;
base_parameters.grid.pml_size = [20] / base_parameters.grid.resolution_mm;
base_parameters.grid.axisymmetric = 1;
base_parameters.subject_id = subject_id;

% Shortened thermal protocol (see PRESTUS_paramvar_08 / README.md, Methods,
% Analysis 4): 20s stimulation, 0s post-stimulation cool-off.
base_parameters.timing.ptrd = 20;
base_parameters.timing.post_ptri_dur = 0;

base_parameters.hpc.timelimit = '00:45:00';  % fresh evaluations took up to 29 min (PRESTUS_paramvar_08);
                                              % 00:15:00 and 00:30:00 both had jobs killed by SLURM

%% ── Factor lists ─────────────────────────────────────────────────────────

acoustic_factors = build_factor_list(tissue_labels, tissue_suffixes, acoustic_properties, acoustic_suffixes);
thermal_factors  = build_factor_list(tissue_labels, tissue_suffixes, thermal_properties, thermal_suffixes);
all_factors = [acoustic_factors, thermal_factors];
n_acoustic = numel(acoustic_factors);
n_thermal  = numel(thermal_factors);
n_total    = n_acoustic + n_thermal;

%% ── Build the r Morris trajectories (design table + job specs) ─────────────

rng(rng_seed);

n_traj = r;
if pilot_mode
    n_traj = 2;
end

% design_rows has one row per point in every trajectory (not just flip
% events), including step 0 (baseline, tissue/property/direction left
% empty): the follow-up analysis needs every step's output_affix to look
% up its outcome and compute elementary_effect(step) = outcome(step) -
% outcome(step-1) for whichever factor step>0's tissue/property records.
n_rows_total = n_traj * (n_total + 1);
design_rows = table('Size', [n_rows_total 6], ...
    'VariableTypes', {'double', 'double', 'string', 'string', 'string', 'string'}, ...
    'VariableNames', {'trajectory', 'step', 'tissue', 'property', 'direction', 'output_affix'});
row_i = 0;

specs = struct('levels', {}, 'output_affix', {});

for tr = 1:n_traj
    baseline_levels = randi(2, 1, n_total);
    acoustic_perm = randperm(n_acoustic);
    thermal_perm  = n_acoustic + randperm(n_thermal);
    step_order = [acoustic_perm, thermal_perm]; % grouped: acoustic factors flipped first, in random order, then thermal

    levels = baseline_levels;
    for step = 0:n_total
        tissue = ""; property = ""; direction = "";
        if step > 0
            flip_idx = step_order(step);
            levels(flip_idx) = 3 - levels(flip_idx); % flip bottom(1)<->top(2)
            tissue = string(all_factors(flip_idx).tissue);
            property = string(all_factors(flip_idx).property);
            direction = string(ternary(levels(flip_idx) == 2, 'bottom_to_top', 'top_to_bottom'));
        end

        output_affix = combo_affix(all_factors, levels);

        row_i = row_i + 1;
        design_rows(row_i, :) = {tr, step, tissue, property, direction, string(output_affix)};

        if ~any(strcmp({specs.output_affix}, output_affix)) % shared point: submit once
            specs(end+1) = struct('levels', levels, 'output_affix', output_affix); %#ok<AGROW>
        end
    end
end

writetable(design_rows, fullfile(pn.path_sim, 'morris_design.csv'));

fprintf('Morris design: %d trajectories x 15 points = %d evaluations, %d unique jobs%s.\n', ...
    n_traj, n_traj * (n_total + 1), numel(specs), ternary(pilot_mode, ' [PILOT MODE]', ''));

%% ── Submit every evaluation (independent jobs, no caching) ─────────────────

for i = 1:numel(specs)
    parameters = build_parameters(base_parameters, all_factors, params_table, ...
        specs(i).levels, specs(i).output_affix);
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

function parameters = build_parameters(base_parameters, all_factors, params_table, levels, output_affix)
% Full, independent evaluation (fresh acoustic + thermal solve, see header).
    level_names = {'bottom', 'top'};
    parameters = base_parameters;
    for f = 1:numel(all_factors)
        parameters.medium_properties.(all_factors(f).tissue).(all_factors(f).property) = ...
            params_table.(all_factors(f).property)(strcmp(params_table.tissue, all_factors(f).tissue) & ...
                                                     strcmp(params_table.level, level_names{levels(f)}));
    end
    parameters.io.output_affix = output_affix;
end
