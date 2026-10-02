function itrusst_apply_calibration(pn, equipment_name, desired_intensity, config_names)
% Copy a calibrated annular drive (amplitude, phase) into the protocol configs.
%
% Reads the transducer yaml written by the PRESTUS calibration
% (data/calibration/PRESTUS_virtual_parameters/<equipment>-F<focus>mm-I<int>wpercm2.yaml)
% and replaces the 'annular: {...}' block of each data/configs/config_<name>.yaml
% with it. All other lines of the configs are left untouched.
%
% Inputs:
%   pn                - path struct from setup_itrusst_paths (uses pn.configs, pn.data)
%   equipment_name    - calibrated setup, e.g. 'itrusst_protocol'
%   desired_intensity - calibrated intensity [W/cm^2], e.g. 60
%   config_names      - (optional) configs to update; default {equipment_name}

if nargin < 4
    config_names = {equipment_name};
end

cal_dir = fullfile(pn.data, 'calibration', 'PRESTUS_virtual_parameters');
cal_file = dir(fullfile(cal_dir, sprintf('%s-F*mm-I%.0fwpercm2.yaml', equipment_name, desired_intensity)));
if numel(cal_file) ~= 1
    error('Expected one calibration yaml for %s at %g W/cm2 in %s, found %d.', ...
        equipment_name, desired_intensity, cal_dir, numel(cal_file));
end
cal_block = read_annular_block(fullfile(cal_dir, cal_file.name));
cal = yaml.loadFile(fullfile(cal_dir, cal_file.name));

for c = 1:numel(config_names)
    config_file = fullfile(pn.configs, ['config_' config_names{c} '.yaml']);
    lines = readlines(config_file);
    [~, idx] = read_annular_block(config_file);
    lines = [lines(1:idx(1)-1); cal_block; lines(idx(end)+1:end)];
    writelines(lines, config_file);

    % check that the drive now matches the calibration
    cfg = yaml.loadFile(config_file);
    if cfg.transducer.annular.elem_p_amp ~= cal.transducer.annular.elem_p_amp || ...
            cfg.transducer.annular.elem_phase_deg ~= cal.transducer.annular.elem_phase_deg
        error('Calibration was not applied correctly to %s.', config_file);
    end
    fprintf('%s: elem_p_amp = %g, elem_phase_deg = %g (from %s)\n', config_names{c}, ...
        cfg.transducer.annular.elem_p_amp, cfg.transducer.annular.elem_phase_deg, cal_file.name);
end

end

function [block, idx] = read_annular_block(yaml_file)
% the 'annular: {...}' flow mapping directly under 'transducer:' (may span lines)
lines = readlines(yaml_file);
first = find(startsWith(lines, '  annular: {'), 1);
if isempty(first)
    error('No inline ''annular: {...}'' block found in %s.', yaml_file);
end
last = first - 1 + find(contains(lines(first:end), '}'), 1);
idx = first:last;
block = lines(idx);
end
