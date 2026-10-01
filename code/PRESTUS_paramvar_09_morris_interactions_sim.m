clear all; close all; clc;

% PRESTUS_PARAMVAR_09_MORRIS_INTERACTIONS_SIM
%
% Benchmark 2 (brain+skull) ONLY. A Morris elementary-effects (one-factor-
% at-a-time trajectory) screening design across all 14 tissue properties,
% specifically to measure INTERACTION strength per property - the gap
% PRESTUS_paramvar_08_constellation_factorial_sim.m's marginal (acoustic-
% exhaustive x thermal-exhaustive-at-one-anchor) design leaves open. This
% is the approach flagged earlier (alongside fractional-factorial and LHS+
% surrogate designs) as the standard, budget-efficient tool for this
% question, revisited here at the user's request rather than a bespoke
% extension of the Analysis-4 design.
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
% ACOUSTIC/THERMAL CACHING (same mechanism as PRESTUS_paramvar_08; see
% doc/PRESTUS_itrusst_methods.md, Analysis 4, for the underlying
% rationale): within each trajectory, the 6 acoustic-relevant factors
% (sound_speed/density/alpha_coeff x brain/skull) are flipped FIRST (in
% random relative order among themselves), then the 8 thermal-relevant
% factors (thermal_conductivity/specific_heat_capacity/perfusion/
% absorption_fraction x brain/skull) - a deliberate grouped ordering (a
% valid Morris design choice; nothing requires a fully unconstrained
% random factor order). This means only the first 7 points of each
% trajectory (baseline + 6 acoustic flips) need a fresh k-Wave acoustic
% solve; the remaining 8 (all-thermal-flips) points reuse the acoustic
% cache from point 7, needing only the cheaper thermal-only re-solve. Every
% point (all 15) runs the thermal stage, since elementary effects on
% riseT_brain/riseT_skull require heating output throughout, not just
% Isppa_brain/Ipa_target.
%
% JOB SUBMISSION: every evaluation is its own job, submitted via
% prestus_pipeline_start(parameters) with parameters.platform set to
% whatever hpc_detect_system() finds ('slurm'/'qsub' if a scheduler is on
% PATH, else falls back to running locally/synchronously) - "auto mode",
% the same dispatcher every other script in this project uses, no custom
% batching or external submission wrapper. All 7 acoustic-cost evaluations
% per trajectory are submitted first (phase 1: 7*r jobs); if running on an
% actual scheduler, this script then BLOCKS (via the toolbox's own
% hpc_wait_for_completion, polling squeue/sacct or qstat every 20s) until
% every phase-1 job has finished, before submitting the 8*r thermal-only
% evaluations (phase 2), each of which depends on its own trajectory's
% acoustic cache already being written. On 'matlab' (no scheduler
% detected), each prestus_pipeline_start call already runs synchronously,
% so phase 1 is inherently complete before phase 2's loop starts.
%
% Total evaluations per trajectory: 15 (7 acoustic-cost + 8 thermal-cost).
% At r=10 (default): 150 jobs total (70 phase-1 + 80 phase-2), ~2.7h of
% compute (see PRESTUS_paramvar_08's header for the timing reference this
% is grounded in), comfortably under the study's <=1000 job-count budget.
% Increase r for more reliable mu*/sigma estimates (cost scales linearly,
% not exponentially - this is the whole point of a screening design over
% a full factorial).
%
% NOT deduplicated across trajectories: with r>~8, different trajectories'
% acoustic-cost points will occasionally land on the same 6-bit acoustic
% combo by chance (pigeonhole: e.g. 70 draws from only 64 possible combos
% at r=10). Every acoustic-cost evaluation always requests its own fresh
% solve regardless (parameters.modules.run_acoustic_sims=1); PRESTUS's own
% io.overwrite_files='never' cache-file check opportunistically skips the
% solve and loads an equivalent, already-written cache instead when
% another evaluation happens to have completed first - a free efficiency
% bonus with no correctness cost either way, since no evaluation ever
% depends on that reuse actually happening.
%
% STORAGE: same reductions as PRESTUS_paramvar_08 (modules.
% run_nifti_creation / generate_report disabled; io.save_thermal_matrices
% disabled; io.save_acoustic_matrices enabled only for acoustic-cost
% evaluations).
%
% ANALYSIS (not yet written - no results exist to analyse): once this
% design has been run, mu*/sigma per factor can be computed directly from
% the trajectory design table this script writes to figures/paramvar/
% morris_design.csv (columns: trajectory, step, tissue, property,
% direction, output_affix) joined against each output_affix's riseT_
% brain/riseT_skull/Isppa_brain/Ipa_target from its sub-002_phantom
% <affix>.csv.
%
% ALWAYS validate with pilot_mode = true (restricts to 2 trajectories)
% before committing the full design - see PRESTUS_paramvar_08's header for
% the same non-decoupled/non-cached cross-check this design needs too.
%
% Output: data/tussim/<variant_name>/sub-002/sub-002_phantom<affix>.csv
% per evaluation; figures/paramvar/morris_design.csv (trajectory design
% table for the follow-up analysis).

user = 'julkos';

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
rootpath = fileparts(pathstr);

addpath(fullfile(rootpath, 'code'));
pn = setup_itrusst_paths(rootpath, user);

%% ── Configuration ───────────────────────────────────────────────────────────

subject_id   = 2;      % Benchmark 2 (brain+skull) only
setup_name   = 'itrusst_protocol1';
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
base_parameters.io.overwrite_files = 'never';   % lets an already-cached acoustic combo be loaded for free if another
                                                 % evaluation happens to finish first (opportunistic only, see header NOTE)
base_parameters.io.overwrite_simnibs = 0;
base_parameters.io.save_heatingvideo = 0;
base_parameters.io.save_matrices = 0;
base_parameters.io.save_thermal_matrices = 0;
base_parameters.modules.run_nifti_creation = 0;
base_parameters.modules.generate_report = 0;
base_parameters.grid.resolution_mm = 0.5;
base_parameters.grid.default_dims = [70, 120] / base_parameters.grid.resolution_mm;
base_parameters.grid.pml_size = [20] / base_parameters.grid.resolution_mm;
base_parameters.grid.axisymmetric = 1;
base_parameters.subject_id = subject_id;

% Shortened thermal protocol (see PRESTUS_paramvar_08 / doc/PRESTUS_itrusst_
% methods.md, Analysis 4): 20s stimulation, 0s post-stimulation cool-off.
base_parameters.timing.ptrd = 20;
base_parameters.timing.post_ptri_dur = 0;

base_parameters.hpc.timelimit = '00:30:00';  % NOT a "generous margin" at 00:15:00 as originally set:
                                              % this pilot run had 4/14 completed evaluations killed by
                                              % SLURM at 80-99% through the thermal stage (see PRESTUS_
                                              % paramvar_10_morris_analysis_plot.m's "partial" report)

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

specs_phase1 = struct('levels', {}, 'acoustic_cache_affix', {}, 'output_affix', {});
specs_phase2 = struct('levels', {}, 'acoustic_cache_affix', {}, 'output_affix', {});

for tr = 1:n_traj
    baseline_levels = randi(2, 1, n_total);
    acoustic_perm = randperm(n_acoustic);
    thermal_perm  = n_acoustic + randperm(n_thermal);
    step_order = [acoustic_perm, thermal_perm]; % grouped: acoustic factors flipped first, in random order, then thermal

    levels = baseline_levels;
    acoustic_affix_at_point7 = ''; % set once all 6 acoustic factors have been flipped
    for step = 0:n_total
        tissue = ""; property = ""; direction = "";
        if step > 0
            flip_idx = step_order(step);
            levels(flip_idx) = 3 - levels(flip_idx); % flip bottom(1)<->top(2)
            tissue = string(all_factors(flip_idx).tissue);
            property = string(all_factors(flip_idx).property);
            direction = string(ternary(levels(flip_idx) == 2, 'bottom_to_top', 'top_to_bottom'));
        end

        acoustic_affix = combo_affix(acoustic_factors, levels(1:n_acoustic));
        output_affix   = combo_affix(all_factors, levels);

        row_i = row_i + 1;
        design_rows(row_i, :) = {tr, step, tissue, property, direction, string(output_affix)};

        if step <= n_acoustic
            specs_phase1(end+1) = struct('levels', levels, 'acoustic_cache_affix', acoustic_affix, ...
                'output_affix', output_affix); %#ok<AGROW>
            if step == n_acoustic
                acoustic_affix_at_point7 = acoustic_affix;
            end
        else
            specs_phase2(end+1) = struct('levels', levels, 'acoustic_cache_affix', acoustic_affix_at_point7, ...
                'output_affix', output_affix); %#ok<AGROW>
        end
    end
end

writetable(design_rows, fullfile(pn.figures, 'morris_design.csv'));

fprintf('Morris design: %d trajectories x 15 points = %d evaluations (%d acoustic-cost, %d thermal-cost)%s.\n', ...
    n_traj, n_traj * (n_total + 1), numel(specs_phase1), numel(specs_phase2), ternary(pilot_mode, ' [PILOT MODE]', ''));

%% ── Phase 1: submit every acoustic-cost evaluation, then wait ──────────────

phase1_job_ids = {}; % cell array: job_id is numeric for SLURM but a string for qsub
for i = 1:numel(specs_phase1)
    parameters = build_parameters(base_parameters, all_factors, params_table, ...
        specs_phase1(i).levels, true, specs_phase1(i).acoustic_cache_affix, specs_phase1(i).output_affix);
    job_id = prestus_pipeline_start(parameters);
    if ~isempty(job_id)
        phase1_job_ids{end+1} = job_id; %#ok<AGROW>
    end
end

if ~strcmp(hpc_type, 'matlab') && ~isempty(phase1_job_ids)
    fprintf('Waiting for %d phase-1 jobs to complete before submitting phase 2...\n', numel(phase1_job_ids));
    for j = 1:numel(phase1_job_ids)
        hpc_wait_for_completion(phase1_job_ids{j}, hpc_type);
    end
end

%% ── Phase 2: submit every thermal-only evaluation (all caches now exist) ───

for i = 1:numel(specs_phase2)
    parameters = build_parameters(base_parameters, all_factors, params_table, ...
        specs_phase2(i).levels, false, specs_phase2(i).acoustic_cache_affix, specs_phase2(i).output_affix);
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

function parameters = build_parameters(base_parameters, all_factors, params_table, ...
    levels, run_fresh_acoustic, acoustic_cache_affix, output_affix)
    level_names = {'bottom', 'top'};
    parameters = base_parameters;
    for f = 1:numel(all_factors)
        parameters.medium_properties.(all_factors(f).tissue).(all_factors(f).property) = ...
            params_table.(all_factors(f).property)(strcmp(params_table.tissue, all_factors(f).tissue) & ...
                                                     strcmp(params_table.level, level_names{levels(f)}));
    end

    parameters.modules.run_source_setup  = run_fresh_acoustic;
    parameters.modules.run_acoustic_sims = run_fresh_acoustic;
    parameters.modules.run_heating_sims  = 1;   % every point needs heating output
    parameters.io.acoustic_cache_affix   = acoustic_cache_affix;
    parameters.io.output_affix           = output_affix;
    parameters.io.save_acoustic_matrices = run_fresh_acoustic;
end
