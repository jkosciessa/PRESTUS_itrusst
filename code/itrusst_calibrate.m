function [opt_source_amp, opt_source_phase_deg, opt_source_phase_rad] = ...
    itrusst_calibrate(pn, parameters, subject_id, desired_intensity, transducer_name)

    % In this function, we run an initial free-water simulation to generate
    % a reference profile that we want to optimize the transducer
    % calibration for.

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
    xlabel('Axial Position [mm]');
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
    xline(opt_res_params.transducer(1).focal_distance_ep, '--');
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

    correctionFactor = (desired_intensity./max(axial_pressure));
    profile_empirical.axial_distance_bowl = pos_x_sim_res;
    profile_empirical.axial_intensity = axial_pressure.*correctionFactor;

    cd(pn.configs);

    parameters.calibration = yaml.loadFile('calibration_config.yaml');
    parameters.calibration.opt_phase_precession = 'monotonic';
    parameters.calibration.opt_amp_validation = 'always';
    parameters.calibration.path_output = fullfile(pn.configs,'..','calibration'); 
    mkdir(parameters.calibration.path_output)
    parameters.calibration.path_output_profiles = ...
        fullfile(parameters.calibration.path_output,'PRESTUS_virtual_parameters');
        mkdir(parameters.calibration.path_output_profiles)

    job_id = calibration_pipeline_start(profile_empirical,...
        transducer_name, ...
        desired_intensity, ...
        parameters.transducer(1).focal_distance_ep, ...
        parameters, ...
        subject_id);

end