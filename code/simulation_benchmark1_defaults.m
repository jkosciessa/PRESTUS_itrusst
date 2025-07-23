% This script executes PRESTUS stimulations for multiple config files
% differing in their parameter assumptions.

clear all; close all; clc;

user = 'julkos' % julkos/marwim/sebrei

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code'));

pn.tuSIM = fullfile(rootpath, 'tools', 'PRESTUS'); addpath(pn.tuSIM);
pn.tuSIM_fun = fullfile(pn.tuSIM, 'functions'); addpath(pn.tuSIM_fun);
pn.tuSIM_tools = fullfile(pn.tuSIM, 'toolboxes'); addpath(genpath(pn.tuSIM_tools));
pn.kwave = fullfile(pn.tuSIM_tools, 'k-wave', 'k-Wave'); addpath(pn.kwave);
if strcmp(user, 'julkos')
    pn.simnibs = fullfile('/home', 'neuromod', 'julkos', '.conda', 'envs', 'simnibs_env', 'bin');
elseif strcmp(user, 'marwim')
    pn.simnibs = fullfile('/home', 'neuromod', 'marwim', 'miniconda3', 'envs', 'simnibs_env', 'bin');
elseif strcmp(user, 'sebrei')
    pn.simnibs = fullfile('/home', 'neuromod', 'sebrei', '.conda', 'envs', 'simnibs_env', 'bin');
end
pn.minimize = fullfile(pn.tuSIM_tools, 'FEX-minimize'); addpath(pn.minimize);
pn.configs = fullfile(rootpath, 'data', 'configs');
pn.data_path = fullfile(rootpath, 'data', 'bids');
pn.data_seg = fullfile(rootpath, 'data', 'simnibs');
pn.nifti = (fullfile(rootpath, 'tools', 'nifti_toolbox')); addpath(pn.nifti);

%% define variables to loop

% setup_list = {
%     'itrusst_protocol1_bottom'; ...
%     'itrusst_protocol1_top'; ...
%     'itrusst_protocol1_mid'; ...
%     'itrusst_protocol1_PRESTUS';...
% 'itrusst_protocol1_kPlan'; 'itrusst_protocol1_BabelBrain'};
setup_list = {'itrusst_protocol1_BabelBrain'};

all_subjects = [001:002];
intensities = 60;
copy_positions = 1;

for subject_id = all_subjects

    for i_setup = 1:length(setup_list)
    for i_intensity = 1:length(intensities)
        transducer_name = setup_list{i_setup};
        desired_intensity = intensities(i_intensity);

        %% load parameter file and adjust paths if necessary
        
        cd(fullfile(rootpath, 'data')) % needs to contain a folder called "configs" in this setup
        parameters = load_parameters(['config_',transducer_name,'.yaml']);

        if ismac
            % reset paths to local mac install; NOTE: horrible not to have this be dynamic
            parameters.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
        else
            parameters.simnibs_bin_path = pn.simnibs;
        end
        pn.data_sims = fullfile(rootpath, 'data', 'tussim', [transducer_name]);
            if ~exist(pn.data_sims); mkdir(pn.data_sims); end
        parameters.ld_library_path ="/opt/gcc/7.2.0/lib64";
        parameters.data_path = pn.data_seg; % use simnibs folder
        parameters.seg_path = pn.data_seg;
        parameters.t1_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
        parameters.t2_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
        parameters.sim_path = pn.data_sims;
        parameters.paths_to_add = {pn.kwave, pn.minimize};
        pn.outputs_folder = fullfile(parameters.sim_path,sprintf('sub-%03d', subject_id));
        if ~exist(pn.outputs_folder); mkdir(pn.outputs_folder); end
        parameters.usepseudoCT = 0;

        %% run free-water simulations

        parameters.simulation_medium = 'water'; % indicate that we only want the simulation in the water medium for now
        parameters.interactive = 0;
        parameters.overwrite_files = 'always';
        parameters.run_heating_sims = 0;
        if ismac % if run locally on a mac
            parameters.using_donders_hpc = 0;
            parameters.code_type = 'matlab_cpu';
        else
            parameters.using_donders_hpc = 1;
        end
        
        single_subject_pipeline(subject_id, parameters);

        %% plot free-water results

        opt_res = load(fullfile(pn.outputs_folder, sprintf('sub-%03d_water_results%s.mat',...
        subject_id, parameters.results_filename_affix)),'sensor_data','parameters');

        p_max = gather(opt_res.sensor_data.p_max_all);
        pred_axial_pressure_opt = squeeze(p_max(opt_res.parameters.transducer.pos_grid(1),:));
        axial_pressure = pred_axial_pressure_opt.^2/...
            (2*parameters.medium.water.sound_speed*parameters.medium.water.density).* 1e-4;

        h = figure('Position', [10 10 900 500]);
        hold on
        xlabel('Axial Position [mm]');
        ylabel('Intensity [W/cm^2]');
        pos_x_axis = (1:parameters.default_grid_dims(2)).*parameters.grid_step_mm; % x-axis [mm]
        pos_x_trans = (opt_res.parameters.transducer.pos_grid(2)-1)*parameters.grid_step_mm; % x-axis position of transducer [mm]
        pos_x_sim_res = pos_x_axis-pos_x_trans; % axial position for the simulated results, relative to transducer position [mm]
        plot(pos_x_sim_res, axial_pressure);
        hold off
        xline(opt_res.parameters.expected_focal_distance_mm, '--');
        yline(desired_intensity, '--');
        plotname = fullfile(pn.outputs_folder, 'simulation_analytic.png');
        saveas(h, plotname, 'png');
        close(h);

        %% introduce a simple free-water amplitude scaling
        % Assumption: the "real" profile linearly scales to a peak of the desired intensity

        % axial_pressure | profile measured in simulation grid
        %                pos_x_axis - raw simulation axis, 1: transducer bowl
        %                pos_x_sim_res - axis without transducer position

        parameters.correctEPdistance = 0; % do not correct for exit plan distance

        correctionFactor = (desired_intensity./max(axial_pressure));
        real_profile(:,1) = pos_x_sim_res; %pos_x_axis;
        real_profile(:,2) = axial_pressure.*correctionFactor;
%         real_profile(real_profile(:,1)<=0, :) = [];

        [opt_source_amp, opt_phases] = transducer_calibration(...
            pn, ...
            parameters, ...
            transducer_name, ...
            subject_id, ...
            desired_intensity, ...
            real_profile);

        %% run soft-tissue simulation

        parameters = load_parameters(['config_',transducer_name,'.yaml']);

        parameters.ld_library_path ="/opt/gcc/7.2.0/lib64";
        parameters.data_path = pn.data_seg;
        parameters.seg_path = pn.data_seg;
        parameters.t1_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
        parameters.t2_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
        parameters.sim_path = pn.data_sims;
        parameters.paths_to_add = {pn.kwave, pn.minimize};

        parameters.simulation_medium = 'phantom'; % use default grid in combo with (final) tissue mask
        parameters.interactive = 0;
        parameters.overwrite_files = 'always';
        parameters.overwrite_simnibs = 0;
        parameters.heatingvideo = 0;

        parameters.transducer.source_amp = opt_source_amp; % use calibrated input intensity

        parameters.transducer.pos_t1_grid = [3, 35];
        parameters.focus_pos_t1_grid = [64, 35];

        parameters.run_source_setup = 1;
        parameters.run_acoustic_sims = 1;
        parameters.run_heating_sims = 1;

        parameters.thermal.n_trials = 400; % 80 s
        parameters.thermal.duty_cycle = 0.1;
        parameters.thermal.stim_duration = 0.2; 
        parameters.thermal.sim_time_steps = 0.02;
        parameters.thermal.post_stim_dur = 400;
        parameters.thermal.post_time_steps = 1;
        parameters.thermal.pri_duration = 0.2;
        parameters.thermal.equal_steps = 1;
        parameters.thermal.cem43_iso = 1;

        parameters.code_type = 'matlab_gpu';

        if ismac % if run locally on a mac
            parameters.using_donders_hpc = 0;
        else
            parameters.using_donders_hpc = 1;
        end
        
        % Run the simulation
        single_subject_pipeline_with_slurm(subject_id, parameters, 0, '00:30:00', 8);

    end
    end
end
