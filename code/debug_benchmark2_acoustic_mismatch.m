% DEBUG_BENCHMARK2_ACOUSTIC_MISMATCH
%
% Probes potential sources of baseline acoustic mismatch between PRESTUS
% (kSpaceFirstOrderAS, axisymmetric 2D) and a full-3D viscoelastic FDTD
% solver (e.g. BabelBrain) for benchmark 2 (skull+brain phantom, sub-002).
%
% Both platforms use the same mid-range tissue properties, so the mismatch
% must arise from one or more of:
%   A. Grid resolution / numerical dispersion
%   B. PML size / boundary reflections
%   C. Transducer source model (kWaveArray discretization)
%   D. Shear-wave physics in skull (present in FDTD, absent in k-wave)
%   E. Axisymmetric vs. 3D geometry
%
% APPROACH (PRESTUS-only, no data from the other platform needed):
%   1. Load the existing baseline result.
%   2. Run PRESTUS at a finer grid resolution to bound dispersion error.
%   3. Run PRESTUS with a larger PML to bound boundary effects.
%   4. Compute an analytical transfer-matrix transmission estimate for the
%      skull layer, then compare the fraction it can explain vs. the
%      observed gap.  Optionally estimate the shear-wave correction.
%   5. Print a tidy summary table.
%
% NOTE: Parts 2-3 submit new simulations. Set RUN_SIMS=false to skip them
%       and only run the analytical part on existing results.

clear all; close all; clc;

RUN_SIMS = false;  % set true to submit new PRESTUS runs
user     = 'julkos';

% -------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'))
rootpath = pwd;

% On macOS the volume may be mounted under /Volumes rather than /project.
% Override rootpath if the derived path does not exist on disk.
if ~isfolder(fullfile(rootpath, 'data'))
    alt = strrep(rootpath, '/project/', '/Volumes/');
    if isfolder(fullfile(alt, 'data'))
        rootpath = alt;
        cd(rootpath);
    end
end

addpath(fullfile(rootpath, 'code'));
pn = setup_itrusst_paths(rootpath, user);

subject_id = 2;
sub_str    = sprintf('sub-%03d', subject_id);

baseline_dir = fullfile(rootpath, 'data', 'tussim', 'medium_mid', sub_str);

% -------------------------------------------------------------------------
% 1.  Load baseline scalar result
% -------------------------------------------------------------------------
baseline_csv = fullfile(baseline_dir, [sub_str '_phantom_output_table.csv']);
if ~isfile(baseline_csv)
    error('Baseline CSV not found:\n  %s\nRun the benchmark sim first.', baseline_csv);
end

T0 = readtable(baseline_csv);
Isppa_brain_baseline = T0.Isppa_brain(1);
focal_z_vox          = T0.focus_pos_vox_3(1);
tx_x                 = T0.trans_pos_vox_1(1);
tx_y                 = T0.trans_pos_vox_2(1);

fprintf('=== BASELINE (itrusst_protocol1a) ===\n');
fprintf('  Isppa_brain          = %.4f W/cm2\n',  Isppa_brain_baseline);
fprintf('  real_focal_dist_mm   = %.2f mm\n',      T0.real_focal_distance_mm(1));
fprintf('  Isppa_skull          = %.4f W/cm2\n',  T0.Isppa_skull(1));

% -------------------------------------------------------------------------
% 2.  Numerical sensitivity sweeps (optional re-runs)
% -------------------------------------------------------------------------
if RUN_SIMS

    % Base parameters
    cd(fullfile(rootpath, 'data', 'configs'));
    parameters_base = load_parameters('config_itrusst_protocol1.yaml');
    parameters_base.path.anat             = pn.path_seg;
    parameters_base.path.seg              = pn.path_seg;
    parameters_base.path.t1_pattern       = fullfile(sprintf('m2m_sub-%03d', subject_id), 'T1.nii.gz');
    parameters_base.path.t2_pattern       = fullfile(sprintf('m2m_sub-%03d', subject_id), 'T2_reg.nii.gz');
    parameters_base.startup.paths_to_add  = {pn.kwave, pn.minimize, ...
        fullfile(pn.prestus_fun, 'telemetry')};
    parameters_base.simulation.medium     = 'phantom';
    parameters_base.simulation.interactive = 0;
    parameters_base.io.overwrite_files    = 'always';
    parameters_base.io.overwrite_simnibs  = 0;
    parameters_base.io.save_heatingvideo  = 0;
    parameters_base.io.save_matrices      = 0;
    parameters_base.simulation.code_type  = 'cpp_cpu';
    parameters_base.simulation.precision  = 'single';
    parameters_base.subject_id            = subject_id;

    if ismac
        parameters_base.startup.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
    end

    sim_dir = fullfile(rootpath, 'data', 'tussim', 'diagnostic_benchmark2');
    if ~exist(sim_dir, 'dir'); mkdir(sim_dir); end

    % Variant A — finer grid resolution (numerical dispersion)
    % default_dims must scale inversely with resolution to keep the same
    % physical domain size (baseline: [140 240] @ 1 mm = 140x240 mm).
    params_res = parameters_base;
    params_res.path.sim        = sim_dir;
    params_res.io.output_affix = '_diag_res050';  % 0.5 mm — less extreme step
    params_res.grid.resolution_mm = 0.5;
    params_res.grid.default_dims  = [280 480];    % same physical domain: 140x240 mm
    fprintf('\n--- Submitting variant: diag_res050 (resolution 0.5 mm) ---\n');
    prestus_pipeline_start(params_res);

    % Variant B — larger PML (boundary reflections)
    params_pml = parameters_base;
    params_pml.path.sim        = sim_dir;
    params_pml.io.output_affix = '_diag_pml60';
    params_pml.grid.pml_size   = 60;
    fprintf('\n--- Submitting variant: diag_pml60 (PML 60 pts) ---\n');
    prestus_pipeline_start(params_pml);
end

% -------------------------------------------------------------------------
% 3.  Load all available diagnostic results and compare to baseline
% -------------------------------------------------------------------------
fprintf('\n=== NUMERICAL SENSITIVITY RESULTS ===\n');

diag_dir = fullfile(rootpath, 'data', 'tussim', 'diagnostic_benchmark2', sub_str);

affixes      = {'diag_res050', 'diag_pml60'};
descriptions = {'resolution 0.5 mm (grid dispersion)', 'PML 60 pts (boundary reflections)'};

fprintf('%-40s  %8s  %8s  %8s\n', 'Variant', 'Isppa_br', 'Delta', 'Delta%');
fprintf('%-40s  %8.4f\n', 'Baseline (1 mm, PML 40)', Isppa_brain_baseline);

for i = 1:numel(affixes)
    fn = fullfile(diag_dir, [sub_str '_phantom_' affixes{i} '.csv']);
    if ~isfile(fn)
        fprintf('  %-38s  (not yet run)\n', descriptions{i});
        continue;
    end
    T = readtable(fn);
    v = T.Isppa_brain(1);
    fprintf('%-40s  %8.4f  %+8.4f  %+7.1f%%\n', ...
        descriptions{i}, v, v - Isppa_brain_baseline, ...
        100*(v - Isppa_brain_baseline)/Isppa_brain_baseline);
end

% -------------------------------------------------------------------------
% 4.  Analytical transfer-matrix: skull transmission + shear estimate
% -------------------------------------------------------------------------
fprintf('\n=== ANALYTICAL TRANSMISSION ESTIMATE ===\n');

% --- tissue properties (mid values from config) ---
c_w   = 1500;  rho_w  = 994;   % water/coupling
c_b   = 1546;  rho_b  = 1045;  % brain (soft tissue)
c_sk  = 2800;  rho_sk = 1850;  % skull (cortical bone, compressional)
cs_sk = 1500;  % shear speed in cortical bone (literature: ~1500 m/s)
% attenuation: dB/(cm*MHz^n) at 0.5 MHz
alpha_w  = 2.17e-3; n_w  = 2;   % water
alpha_b  = 0.42;    n_b  = 1;   % brain
alpha_sk = 13.3;    n_sk = 1;   % skull
freq_MHz = 0.5;

% --- skull layer thickness from medium mask ---
mask_fn = fullfile(baseline_dir, [sub_str '_final_medium_masks_orig_coord.nii.gz']);

% Read actual simulation resolution from the log rather than the NIfTI
% header: prior to the nifti_to_t1w fix, phantom outputs were always
% written with pixdim=1 regardless of grid.resolution_mm.
log_files = dir(fullfile(baseline_dir, 'batch_job_logs', '*.log'));
sim_res_mm = NaN;
for li = 1:numel(log_files)
    txt = fileread(fullfile(log_files(li).folder, log_files(li).name));
    tok = regexp(txt, 'resolution_mm:\s+([\d.]+)', 'tokens', 'once');
    if ~isempty(tok)
        sim_res_mm = str2double(tok{1});
        break;
    end
end
if isnan(sim_res_mm)
    warning('Could not read resolution from log; assuming 0.5 mm.');
    sim_res_mm = 0.5;
end
fprintf('  Simulation resolution (from log): %.2f mm\n', sim_res_mm);

skull_thickness_cm = NaN;
if isfile(mask_fn)
    mask = niftiread(mask_fn);
    dz_cm = sim_res_mm / 10;  % use config resolution, not NIfTI header

    axial_col = squeeze(mask(tx_x, tx_y, :));
    % skull labels: anything that is not water(1) or brain(2) or skin(3)
    % Identify skull as label with highest sound speed — use label 4 or 5
    % (PRESTUS phantom layers: 1=water, 2=skin, 3=brain, 4=skull_trabecular,
    %  5=skull_cortical, 6=skull_inner; exact mapping depends on config)
    % Here we flag any label >= 4 as bone-like:
    skull_mask_1d = axial_col >= 4;
    n_skull_vox   = sum(skull_mask_1d);
    skull_thickness_cm = n_skull_vox * dz_cm;
    fprintf('  Skull thickness along beam axis: %d vox = %.2f mm\n', ...
        n_skull_vox, skull_thickness_cm * 10);
else
    % fallback: typical parietal bone ~6 mm
    skull_thickness_cm = 0.6;
    warning('Medium mask not found; using assumed skull thickness %.1f mm.', ...
        skull_thickness_cm * 10);
end

% --- plane-wave attenuation in each layer to focal point ---
% Layer thicknesses are approximate; skull thickness from mask above.
% Water/coupling: transducer face to skull surface (~5 mm typical for phantom)
% Skull: as measured above
% Brain: skull surface to focus

% Load focal distance and transducer position from baseline CSV
focal_dist_mm = T0.real_focal_distance_mm(1);      % focus from exit plane
skull_start_mm = 5;                                  % coupling water (estimated)
brain_to_focus_mm = focal_dist_mm - skull_thickness_cm*10 - skull_start_mm;

L_w_cm  = skull_start_mm / 10;
L_sk_cm = skull_thickness_cm;
L_b_cm  = max(brain_to_focus_mm, 0) / 10;

fprintf('  Layer thicknesses used for attenuation estimate:\n');
fprintf('    Water (coupling): %.1f mm\n', L_w_cm*10);
fprintf('    Skull:            %.1f mm\n', L_sk_cm*10);
fprintf('    Brain:            %.1f mm\n', L_b_cm*10);

% Plane-wave attenuation (dB)
dB_w  = alpha_w  * (freq_MHz^n_w)  * L_w_cm;
dB_sk = alpha_sk * (freq_MHz^n_sk) * L_sk_cm;
dB_b  = alpha_b  * (freq_MHz^n_b)  * L_b_cm;
dB_total = dB_w + dB_sk + dB_b;

fprintf('\n  Plane-wave attenuation:\n');
fprintf('    Water: %.4f dB\n',  dB_w);
fprintf('    Skull: %.4f dB\n',  dB_sk);
fprintf('    Brain: %.4f dB\n',  dB_b);
fprintf('    Total: %.4f dB  ->  intensity factor = %.4f\n', ...
    dB_total, 10^(-dB_total/10));

% --- Reflection loss at interfaces (normal incidence) ---
Z_w  = rho_w  * c_w;
Z_b  = rho_b  * c_b;
Z_sk = rho_sk * c_sk;

R_w_sk = ((Z_sk - Z_w) / (Z_sk + Z_w))^2;  % water->skull reflection
T_w_sk = 1 - R_w_sk;
R_sk_b = ((Z_b - Z_sk)  / (Z_b  + Z_sk))^2; % skull->brain reflection
T_sk_b = 1 - R_sk_b;

fprintf('\n  Interface reflection losses (normal incidence):\n');
fprintf('    Water->Skull:  R = %.4f  T_intensity = %.4f\n', R_w_sk, T_w_sk);
fprintf('    Skull->Brain:  R = %.4f  T_intensity = %.4f\n', R_sk_b, T_sk_b);

total_transmission = T_w_sk * T_sk_b * 10^(-dB_total/10);
fprintf('\n  Combined transmission (attenuation + reflection): %.4f\n', total_transmission);

% --- Shear wave estimate (energy converted from compressional to shear in skull) ---
% At normal incidence there is no mode conversion; at oblique incidence
% (relevant for a focused bowl transducer) some energy couples to shear.
% Rule-of-thumb: for 500 kHz and typical parietal bone, mode conversion
% reduces transmitted intensity by an additional ~10-25% relative to
% purely acoustic treatment (from literature on transcranial ultrasound).
fprintf('\n  Shear-wave mode conversion in skull:\n');
fprintf('    Not modelled in kSpaceFirstOrderAS (PRESTUS).\n');
fprintf('    In viscoelastic FDTD (BabelBrain) shear conversion is implicit.\n');
fprintf('    Literature estimate for 500 kHz: ~10-25%% additional attenuation.\n');
fprintf('    -> Expected Isppa_brain reduction: %.1f–%.1f W/cm2 below PRESTUS output.\n', ...
    0.10 * Isppa_brain_baseline, 0.25 * Isppa_brain_baseline);

% -------------------------------------------------------------------------
% 5.  Summary
% -------------------------------------------------------------------------
fprintf('\n=== SUMMARY OF CANDIDATE SOURCES ===\n');
fprintf('%-45s  %s\n', 'Source', 'Testable with PRESTUS alone?');
fprintf('%-45s  %s\n', 'A. Grid dispersion (resolution ~0.5 mm)', 'YES — run diag_res025 variant');
fprintf('%-45s  %s\n', 'B. PML reflections', 'YES — run diag_pml60 variant');
fprintf('%-45s  %s\n', 'C. kWaveArray source discretization', 'Partially — compare kWaveArray on/off');
fprintf('%-45s  %s\n', 'D. Shear-wave mode conversion in skull', 'NO — requires viscoelastic solver');
fprintf('%-45s  %s\n', 'E. Axisymmetric vs 3D geometry', 'YES — switch grid.axisymmetric to 0');
fprintf('\nIf A+B+E are small (< ~2%%), the residual gap is most likely\n');
fprintf('explained by shear-wave physics (D) absent from kSpaceFirstOrderAS.\n');
