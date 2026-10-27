function itrusst_apply_calibration(pn, equipment_name, desired_intensity, config_names)
% Copy a calibrated annular drive (amplitude, phase) into the protocol configs.
%
% Reads the transducer yaml written by the PRESTUS calibration
% (data/calibration/PRESTUS_virtual_parameters/<equipment>-F<focus>mm-I<int>wpercm2.yaml)
% and replaces the 'annular: {...}' block of each data/configs/config_<name>.yaml
% with it. It also sets transducer.target_isppa_wcm2 and transducer.tx_calibration
% to the calibrated intensity, so the simulations model that intensity instead of
% the PRESTUS default target (30 W/cm^2). All other lines are left untouched.
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
if isfield(cal.transducer, 'tx_calibration') && isfield(cal.transducer.tx_calibration, 'isppa_wcm2') && ...
        abs(cal.transducer.tx_calibration.isppa_wcm2 - desired_intensity) > 1e-6
    error('%s records a calibration for %g W/cm2, not %g W/cm2.', cal_file.name, ...
        cal.transducer.tx_calibration.isppa_wcm2, desired_intensity);
end
isppa_lines = [sprintf("  target_isppa_wcm2: %.1f   # free-water ISPPA the annular drive above was calibrated for (itrusst_apply_calibration)", desired_intensity); ...
               sprintf("  tx_calibration: {isppa_wcm2: %.1f, phase_source: library, isppa_source: library}", desired_intensity)];

for c = 1:numel(config_names)
    config_file = fullfile(pn.configs, ['config_' config_names{c} '.yaml']);
    lines = readlines(config_file);
    [~, idx] = read_annular_block(config_file);
    lines = [lines(1:idx(1)-1); cal_block; lines(idx(end)+1:end)];
    % drop earlier target / calibration lines, then add them after the annular block
    lines = lines(~startsWith(lines, ["  target_isppa_wcm2:", "  tx_calibration:"]));
    [~, idx] = read_annular_block_lines(lines, config_file);
    lines = [lines(1:idx(end)); isppa_lines; lines(idx(end)+1:end)];
    writelines(lines, config_file);

    % check that the drive now matches the calibration
    cfg = yaml.loadFile(config_file);
    if cfg.transducer.annular.elem_p_amp ~= cal.transducer.annular.elem_p_amp || ...
            cfg.transducer.annular.elem_phase_deg ~= cal.transducer.annular.elem_phase_deg || ...
            cfg.transducer.target_isppa_wcm2 ~= desired_intensity || ...
            cfg.transducer.tx_calibration.isppa_wcm2 ~= desired_intensity
        error('Calibration was not applied correctly to %s.', config_file);
    end
    fprintf('%s: elem_p_amp = %g, elem_phase_deg = %g, target_isppa_wcm2 = %g (from %s)\n', config_names{c}, ...
        cfg.transducer.annular.elem_p_amp, cfg.transducer.annular.elem_phase_deg, ...
        cfg.transducer.target_isppa_wcm2, cal_file.name);
end

end

function [block, idx] = read_annular_block(yaml_file)
% the 'annular: {...}' flow mapping directly under 'transducer:' (may span lines)
[block, idx] = read_annular_block_lines(readlines(yaml_file), yaml_file);
end

function [block, idx] = read_annular_block_lines(lines, yaml_file)
first = find(startsWith(lines, '  annular: {'), 1);
if isempty(first)
    error('No inline ''annular: {...}'' block found in %s.', yaml_file);
end
last = first - 1 + find(contains(lines(first:end), '}'), 1);
idx = first:last;
block = lines(idx);
end
