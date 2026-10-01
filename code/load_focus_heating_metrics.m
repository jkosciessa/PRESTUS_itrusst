function [results, xlabels] = load_focus_heating_metrics(data_dir, subject_id, ...
    tissue_labels, tissue_suffixes, property_labels, property_suffixes, sweep_values, property_symbols)
% Load heating-rise metrics at two point locations, read directly out of
% the per-condition heatrise NIfTIs, to disambiguate "heating at the
% focus" from the whole-brain independent maximum (riseT_brain, see
% load_paramvar_results.m — its max can land on a different voxel than
% either location below):
%
%   riseT_focus  - temperature rise at the EMPIRICAL acoustic focus voxel
%                  for that specific run (Ix_brain_vox, Iz_brain_vox —
%                  wherever Isppa actually peaked; shifts with tissue
%                  properties)
%   riseT_target - temperature rise at the FIXED nominal target voxel
%                  (focus_pos_vox_1, focus_pos_vox_3 — the config-specified
%                  target; the same voxel for every condition)
%
% Assumes 2-D axisymmetric phantoms [x (lateral), z (axial)] as produced by
% PRESTUS_paramvar_01_benchmarks_sim.m (no singleton y-dimension — see
% PRESTUS_paramvar_04_axial_profiles_2D.m), so only the 1st/3rd voxel
% components of trans_pos_vox/focus_pos_vox apply (the 2nd is a constant
% placeholder equal to the lateral centre, unused here).
%
% Inputs / outputs follow the same convention as load_paramvar_results.m.

if nargin < 8
    property_symbols = [];
end

n_conds = length(tissue_labels) * length(property_labels);
n_sweep = length(sweep_values);
sub_str = sprintf('sub-%03d', subject_id);
nii_dir = fullfile(data_dir, 'nii');

results = struct('riseT_focus', nan(n_sweep, n_conds), 'riseT_target', nan(n_sweep, n_conds));
xlabels = cell(1, n_conds);

required_cols = {'Ix_brain_vox', 'Iz_brain_vox', 'focus_pos_vox_1', 'focus_pos_vox_3'};

i_cond = 0;
for t = 1:length(tissue_labels)
    tissue_suffix = tissue_suffixes{t};

    for p = 1:length(property_labels)
        i_cond = i_cond + 1;
        prop_suffix = property_suffixes{p};
        tissue_sup = tissue_suffix(2:end);
        if ~isempty(property_symbols)
            xlabels{i_cond} = [property_symbols{p} '_{' tissue_sup '}'];
        else
            xlabels{i_cond} = [tissue_sup '_' property_labels{p}];
        end

        for v = 1:length(sweep_values)
            affix    = [tissue_suffix prop_suffix num2str(sweep_values(v))];
            csv_file = fullfile(data_dir, [sub_str '_phantom' affix '.csv']);
            nii_file = fullfile(nii_dir, [sub_str '_phantom_T1w' affix '_heatrise.nii.gz']);

            if ~isfile(csv_file)
                warning('File not found: %s', csv_file);
                continue;
            end
            if ~isfile(nii_file)
                warning('Heatrise NIfTI not found: %s', nii_file);
                continue;
            end

            T = readtable(csv_file);
            if ~all(ismember(required_cols, T.Properties.VariableNames))
                warning('File %s missing focus/target voxel columns.', csv_file);
                continue;
            end

            heatrise = double(niftiread(nii_file));

            ix_focus  = T.Ix_brain_vox(1);
            iz_focus  = T.Iz_brain_vox(1);
            ix_target = T.focus_pos_vox_1(1);
            iz_target = T.focus_pos_vox_3(1);

            results.riseT_focus(v, i_cond)  = heatrise(ix_focus, iz_focus);
            results.riseT_target(v, i_cond) = heatrise(ix_target, iz_target);
        end
    end
end

end
