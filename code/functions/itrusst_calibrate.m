function [opt_source_amp, opt_source_phase_deg, opt_source_phase_rad] = ...
    itrusst_calibrate(pn, parameters, subject_id, desired_intensity, transducer_name)

    % Run a free-water simulation with the configured drive, then rescale
    % elem_p_amp so that the simulated on-axis peak intensity equals
    % desired_intensity (single-element transducers only). Writes
    % data/calibration/PRESTUS_virtual_parameters/<name>-F<focus>mm-I<int>wpercm2.yaml,
    % which itrusst_apply_calibration copies into the protocol config.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run free-water simulations %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    parameters.simulation.medium = 'water'; % indicate that we only want the simulation in the water medium for now
    parameters.simulation.interactive = 0;
    parameters.simulation.code_type = 'cpp_cpu';
    parameters.io.overwrite_files = 'always';
    parameters.io.save_acoustic_matrices = 1; % we want to load the output below
    parameters.hpc.wait_for_job = 1;
    parameters.modules.run_heating_sims = 0;
    % simulate the configured drive as it is: with a target ISPPA and no matching
    % tx_calibration, PRESTUS would scale the field to the target (APPLY_ISPPA_SCALING)
    % and the calibration below would always find the target intensity
    parameters.transducer(1).target_isppa_wcm2 = [];
    if isfield(parameters.transducer, 'tx_calibration')
        parameters.transducer = rmfield(parameters.transducer, 'tx_calibration');
    end
    parameters.modules.run_water_baseline = 0;
    if ~isfield(parameters, 'platform')
        parameters.platform = 'auto';
    end

    parameters.subject_id = subject_id;
    prestus_pipeline_start(parameters);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot free-water results %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cache_dir = fullfile(pn.outputs_folder, 'cache');
    opt_res = load(fullfile(cache_dir, sprintf('sub-%03d_water_results%s.mat',...
        subject_id, parameters.io.output_affix)));

    opt_res_params = opt_res.acoustic_info.parameters;
    if isfield(opt_res.sensor_data, 'isppa_scale_p') && opt_res.sensor_data.isppa_scale_p ~= 1
        error('The free-water field was scaled by %g (ISPPA target); the calibration needs the unscaled field.', ...
            opt_res.sensor_data.isppa_scale_p);
    end

    p_max = gather(opt_res.sensor_data.p_max_all);

    if numel(size(p_max))==2
        pred_axial_pressure_opt = squeeze(p_max(opt_res_params.transducer(1).trans_pos(1),:));
    elseif numel(size(p_max))==3
        pred_axial_pressure_opt = squeeze(p_max(opt_res_params.transducer(1).trans_pos(1),opt_res_params.transducer(1).trans_pos(2), :));
    end
    axial_pressure = pred_axial_pressure_opt.^2/...
        (2*opt_res_params.medium_properties.water.sound_speed*opt_res_params.medium_properties.water.density).* 1e-4;

    h = figure('Position', [10 10 900 500]);
    hold on
    xlabel('Axial distance from bowl apex [mm]');
    ylabel('Intensity [W/cm^2]');
    if numel(size(p_max))==2
        long_dim = 2;
    elseif numel(size(p_max))==3
        long_dim = 3;
    end
    pos_x_axis = (1:size(p_max, long_dim)).*opt_res_params.grid.resolution_mm; % x-axis [mm]
    pos_x_trans = opt_res_params.transducer(1).trans_pos(long_dim)*opt_res_params.grid.resolution_mm; % x-axis position of transducer [mm]
    pos_x_sim_res = pos_x_axis-pos_x_trans; % axial position for the simulated results, relative to transducer position [mm]
    plot(pos_x_sim_res, axial_pressure);
    hold off
    xline(opt_res_params.transducer(1).focal_distance_bowl, '--'); % x-axis is from the bowl apex, not the exit plane
    yline(desired_intensity, '--');
    plotname = fullfile(pn.outputs_folder, 'simulation_analytic.png');
    saveas(h, plotname, 'png');
    saveas(h, strrep(plotname, '.png', '.eps'), 'epsc');
    close(h);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% introduce a simple free-water amplitude scaling %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % The virtual ITRUSST transducer has no equipment entry, so the PRESTUS
    % library calibration (calibration_transducer) does not apply. With a
    % single element the phase has no effect, and in the linear free-water
    % simulation the intensity scales with the square of the source pressure:
    % rescale elem_p_amp so that the simulated on-axis peak equals desired_intensity.

    % scale the drive that was simulated (PRESTUS may adjust the configured one)
    annular = opt_res_params.transducer(1).annular;
    if annular.elem_p_amp ~= parameters.transducer(1).annular.elem_p_amp
        warning('itrusst_calibrate:driveChanged', ...
            'Simulated elem_p_amp (%g Pa) differs from the configured %g Pa; scaling the simulated one.', ...
            annular.elem_p_amp, parameters.transducer(1).annular.elem_p_amp);
    end
    if annular.elem_n ~= 1
        error('itrusst_calibrate rescales the drive of a single-element transducer only; %s has %g elements.', ...
            transducer_name, annular.elem_n);
    end
    [peak_intensity, i_peak] = max(axial_pressure);
    scale = sqrt(desired_intensity ./ peak_intensity);
    fprintf('%s: simulated peak %.2f W/cm2 at %.1f mm from the bowl (geometric focus %.1f mm); elem_p_amp %.1f -> %.1f Pa (x%.5f)\n', ...
        transducer_name, peak_intensity, pos_x_sim_res(i_peak), opt_res_params.transducer(1).focal_distance_bowl, ...
        annular.elem_p_amp, annular.elem_p_amp * scale, scale);

    opt_source_amp       = annular.elem_p_amp * scale;
    opt_source_phase_deg = annular.elem_phase_deg;
    opt_source_phase_rad = deg2rad(opt_source_phase_deg);

    p = parameters;
    p.transducer = parameters.transducer(1);
    p.transducer.annular.elem_p_amp     = opt_source_amp;
    p.transducer.annular.elem_phase_deg = opt_source_phase_deg;
    p.transducer.annular.elem_phase_rad = opt_source_phase_rad;
    p.calibration = struct();
    p.calibration.path_output_profiles      = fullfile(pn.configs, '..', 'calibration', 'PRESTUS_virtual_parameters');
    p.calibration.filename_calibrated_CSV   = ['PRESTUS_virtual_' transducer_name '.csv'];
    p.calibration.equipment_name            = transducer_name;
    p.calibration.desired_focal_distance_ep = parameters.transducer(1).focal_distance_ep;
    p.calibration.desired_intensity         = desired_intensity;
    p.calibration.forward_model             = 'kwave';
    p.calibration.calibration_type          = 'free_water_amplitude_scaling';
    p.calibration.focal_reference           = 'axial_peak'; % else save_optimized_values records flhm_3db
    if ~isfolder(p.calibration.path_output_profiles); mkdir(p.calibration.path_output_profiles); end
    save_optimized_values(p);

    % check the written yaml
    cal_yaml = dir(fullfile(p.calibration.path_output_profiles, ...
        sprintf('%s-F*mm-I%.0fwpercm2.yaml', transducer_name, desired_intensity)));
    if numel(cal_yaml) ~= 1
        error('Expected one calibration yaml for %s at %g W/cm2 in %s, found %d.', ...
            transducer_name, desired_intensity, p.calibration.path_output_profiles, numel(cal_yaml));
    end
    cal = yaml.loadFile(fullfile(cal_yaml.folder, cal_yaml.name));
    if abs(cal.transducer.annular.elem_p_amp - opt_source_amp) > 1e-6 * opt_source_amp
        error('%s holds elem_p_amp = %g, not the computed %g.', cal_yaml.name, ...
            cal.transducer.annular.elem_p_amp, opt_source_amp);
    end

end
