function ref = paper_reference_metrics(data_dir, subject_id)
% Outcomes of the reference medium (all tissue properties at 'mid') for a
% paramvar subject, used as reference lines in the paper range panels.
%
% Reads <sub>_phantom_baseline.csv (PRESTUS_paramvar_00_baseline_sim.m)
% and the empirical-focus heat rise from its heatrise NIfTI (as in
% load_focus_heating_metrics). Without a baseline run, the acoustic metrics
% are taken from a run that varies a thermal property only (brain thermal
% conductivity, bottom level), whose acoustic solution is identical to the
% reference; thermal metrics are then NaN.
%
% Output struct: Ipa_target, Isppa_brain, riseT_skull, riseT_brain,
% riseT_focus, riseT_target (fixed nominal target voxel), has_baseline.

sub_str = sprintf('sub-%03d', subject_id);
ref = struct('Ipa_target', NaN, 'Isppa_brain', NaN, 'riseT_skull', NaN, ...
    'riseT_brain', NaN, 'riseT_focus', NaN, 'riseT_target', NaN, 'has_baseline', false);

base_csv = fullfile(data_dir, [sub_str '_phantom_baseline.csv']);
if isfile(base_csv)
    T = readtable(base_csv);
    ref.has_baseline = true;
    for f = {'Ipa_target', 'Isppa_brain', 'riseT_skull', 'riseT_brain'}
        if ismember(f{1}, T.Properties.VariableNames); ref.(f{1}) = T.(f{1})(1); end
    end
    nii = fullfile(data_dir, 'nii', [sub_str '_phantom_T1w_baseline_heatrise.nii.gz']);
    if isfile(nii) && all(ismember({'Ix_brain_vox', 'Iz_brain_vox'}, T.Properties.VariableNames))
        heatrise = double(niftiread(nii));
        ref.riseT_focus = heatrise(T.Ix_brain_vox(1), T.Iz_brain_vox(1));
        if all(ismember({'focus_pos_vox_1', 'focus_pos_vox_3'}, T.Properties.VariableNames))
            ref.riseT_target = heatrise(T.focus_pos_vox_1(1), T.focus_pos_vox_3(1));
        end
    end
    return;
end

acoustic_csv = fullfile(data_dir, [sub_str '_phantom_b_k_1.csv']);
if isfile(acoustic_csv)
    T = readtable(acoustic_csv);
    ref.Ipa_target  = T.Ipa_target(1);
    ref.Isppa_brain = T.Isppa_brain(1);
    warning('%s: no baseline run, acoustic reference from %s; no thermal reference.', ...
        sub_str, acoustic_csv);
else
    warning('%s: no baseline run and no acoustic substitute in %s.', sub_str, data_dir);
end
end
