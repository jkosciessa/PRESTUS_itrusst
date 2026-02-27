function [opt_source_amp, opt_source_phase_deg, opt_source_phase_rad] = itrusst_calibrate(pn, parameters, subject_id, desired_intensity, transducer_name)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run free-water simulations %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    parameters.simulation_medium = 'water'; % indicate that we only want the simulation in the water medium for now
    parameters.interactive = 0;
    parameters.overwrite_files = 'always';
    parameters.run_heating_sims = 0;
    parameters.code_type = 'matlab_gpu';
    parameters.savemat = 1; % we want to load the output below

    single_subject_pipeline(subject_id, parameters);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot free-water results %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    opt_res = load(fullfile(pn.outputs_folder, sprintf('sub-%03d_water_results%s.mat',...
        subject_id, parameters.results_filename_affix)));

    opt_res_params = opt_res.acoustic_info.parameters;

    p_max = gather(opt_res.sensor_data.p_max_all);

    if numel(size(p_max))==2
        pred_axial_pressure_opt = squeeze(p_max(opt_res_params.transducer.trans_pos(1),:));
    elseif numel(size(p_max))==3
        pred_axial_pressure_opt = squeeze(p_max(opt_res_params.transducer.trans_pos(1),opt_res_params.transducer.trans_pos(2), :));
    end
    axial_pressure = pred_axial_pressure_opt.^2/...
        (2*opt_res_params.medium.water.sound_speed*opt_res_params.medium.water.density).* 1e-4;

    h = figure('Position', [10 10 900 500]);
    hold on
    xlabel('Axial Position [mm]');
    ylabel('Intensity [W/cm^2]');
    if numel(size(p_max))==2
        long_dim = 2;
    elseif numel(size(p_max))==3
        long_dim = 3;
    end
    pos_x_axis = (1:opt_res_params.default_grid_dims(long_dim)).*opt_res_params.grid_step_mm; % x-axis [mm]
    pos_x_trans = (opt_res_params.transducer.trans_pos(long_dim)-1)*opt_res_params.grid_step_mm; % x-axis position of transducer [mm]
    pos_x_sim_res = pos_x_axis-pos_x_trans; % axial position for the simulated results, relative to transducer position [mm]
    plot(pos_x_sim_res, axial_pressure);
    hold off
    xline(opt_res_params.expected_focal_distance_mm, '--');
    yline(desired_intensity, '--');
    plotname = fullfile(pn.outputs_folder, 'simulation_analytic.png');
    saveas(h, plotname, 'png');
    close(h);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% introduce a simple free-water amplitude scaling %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Assumption: the "real" profile linearly scales to a peak of the desired intensity

    % axial_pressure | profile measured in simulation grid
    %                pos_x_axis - raw simulation axis, 1: transducer bowl
    %                pos_x_sim_res - axis without transducer position

    parameters.correctEPdistance = 0; % do not correct for exit plan distance

    correctionFactor = (desired_intensity./max(axial_pressure));
    profile_empirical.dist_from_tran = pos_x_sim_res;
    profile_empirical.profile_focus = axial_pressure.*correctionFactor;
    profile_empirical.focus_wrt_exit_plane = 64; % not really the exit plane here...

    cd(pn.configs);

    parameters.calibration = yaml.loadFile('calibration_config.yaml');
    parameters.calibration.path_output = fullfile(pn.configs,'..','calibration'); 
        mkdir(parameters.calibration.path_output)
    parameters.calibration.path_output_profiles = ...
        fullfile(parameters.calibration.path_output,'PRESTUS_virtual_parameters');
        mkdir(parameters.calibration.path_output_profiles)

    [opt_source_amp, opt_source_phase_deg, opt_source_phase_rad] = ...
        calibration_transducer(profile_empirical,...
        transducer_name, ...
        desired_intensity, ...
        parameters, ...
        subject_id);

end