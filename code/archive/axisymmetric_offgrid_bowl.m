% Simple example of using the axisymmetric code.
%
% author: 
% date: 25th September 2017
% last update: 2nd November 2017

clearvars;
close all

% setup paths

currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;
pn.tuSIM = fullfile(rootpath, 'tools', 'PRESTUS'); addpath(pn.tuSIM);
pn.tuSIM_fun = fullfile(pn.tuSIM, 'functions'); addpath(pn.tuSIM_fun);
pn.tuSIM_tools = fullfile(pn.tuSIM, 'toolboxes'); addpath(genpath(pn.tuSIM_tools));
pn.kwave = fullfile(pn.tuSIM_tools, 'k-wave', 'k-Wave'); addpath(pn.kwave);

% create the computational grid
Nx = 256;           
Ny = 128;          
dx = 0.5e-3;    
dy = dx;          
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% source parameters
SOURCE_DIAMETER = 62.2e-3;  % [m] 
SOURCE_ROC      = 64e-3;  % [m] 
SOURCE_FREQ     = 0.5e6;    % [Hz]
% place the bowl on the edge of the grid so focal point is the first
% lateral grid point
focus_pos_half = [0 , kgrid.y_vec(1)];
focus_pos_full = [0 , 0+eps]; % eps =  smallest step in MATLAB
% place the distributed source outside the PML
x_offset = 30;
y_offset = 0;

% medium parameters
WATER_TEMPERATURE = 20;   % degrees C
DENSITY           = 998; % [kg/m^3]
ALPHA_POWER       = 2;      %
BONA              = 4.986;  %

% compute sound speed and attenuation from temperature
SOUND_SPEED = waterSoundSpeed(WATER_TEMPERATURE);                   % [m/s]
ALPHA_COEFF = waterAbsorption(SOURCE_FREQ/1e6, WATER_TEMPERATURE); % [dB/(MHz^y cm)]

%% Acoustic medium

medium.sound_speed = SOUND_SPEED;
medium.density = DENSITY;
%medium.BonA = BONA;
%medium.alpha_coeff = ALPHA_COEFF;
%medium.alpha_power = ALPHA_POWER;

% create the time array
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed, 0.15);
t_max = 1.2*sqrt(Nx^2 + (2*Ny)^2) * dx / medium.sound_speed;
kgrid.t_array = 0:dt:t_max;

% define time step using an integer number of PPP
PPW = SOUND_SPEED/(SOURCE_FREQ * dx)
%PPP = floor(PPW/CFL);
% dt = 1/(SOURCE_FREQ*PPP);
PPP = 1 / (SOURCE_FREQ * dt)
CFL = PPW / PPP

sensor.mask = ones(Nx, Ny);
sensor.record = {'p_max_all'};

%% create one sided source on edge of grid
% % create empty array
% karray_half = kWaveArray('BLITolerance', 0.02, 'UpsamplingRate', 50);
% % position of the central point on the bowl surface
% position = [kgrid.x_vec(1) + x_offset*dx, kgrid.y_vec(1)];
% % add the arc element (2D bowl)
% karray_half.addArcElement(position, SOURCE_ROC, SOURCE_DIAMETER, focus_pos_half)
% 
% % create a CW signal
% source_signal = createCWSignals(kgrid.t_array, SOURCE_FREQ, 1, 0, 5);
% 
% % make a new kGrid at full size?
% source_half.p = karray_half.getDistributedSourceSignal(kgrid, source_signal);
% 
% source_half.p_mask = karray_half.getArrayBinaryMask(kgrid);
% 
% % run the simulation
% sensor_data_half = kspaceFirstOrderAS(kgrid, medium, source_half, sensor, 'PlotSim', false);
% % 

%% create double size source and crop
% make a new grid with odd number of pts. 
kgrid_mirrored = kWaveGrid(Nx, dx, 2*Ny - 1, dy);

karray_full = kWaveArray('BLITolerance', 0.01, 'UpsamplingRate', 100);
% position of the central point on the bowl surface
position_full = [kgrid.x_vec(1) + x_offset*dx, 0+eps];
% add the arc element (2D bowl)
karray_full.addArcElement(position_full, SOURCE_ROC, SOURCE_DIAMETER, focus_pos_full)

% create a CW signal
source_signal = createCWSignals(kgrid.t_array, SOURCE_FREQ, 1, 0, 5);

%source.p = karray_full.getDistributedSourceSignal(kgrid_mirrored, source_signal);
% get the source mask
source.p_mask = karray_full.getArrayBinaryMask(kgrid_mirrored);
% get the source weighting for each point in the mask
grid_weights = karray_full.getElementGridWeights(kgrid_mirrored, 1);

source.p_mask = source.p_mask(:,Ny:end);
grid_weights = grid_weights(:, Ny:end);

tmp = grid_weights(source.p_mask);
source.p = bsxfun(@times, source_signal, tmp(:));

sensor.mask = ones(Nx, Ny);
sensor.record = {'p_max_all'};

% run the simulation
sensor_data_full = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'PlotSim', false, 'RadialSymmetry', 'WSWA-FFT', 'PMLInside', true);

% plot
figure;
imagesc(sensor_data_full.p_max_all);
axis image

%x = kgrid.x_vec(x_offset:end) - kgrid.x_vec(x_offset);
x = kgrid.x_vec - position_full(1);

figure;
%plot(x,  sensor_data_half.p_max_all(:,1) ./ max(sensor_data_half.p_max_all(:)));
hold all
plot(x,  sensor_data_full.p_max_all(:,1) ./ max(sensor_data_full.p_max_all(:)));

% compute the O'Neil solution and plot for comparison
% positions to compute solution in x propagation direction
x_vec = 0:dx:120e-3;
y_vec = 0:dy:60e-3;
[p_axial, p_lateral] = focusedBowlONeil(SOURCE_ROC, SOURCE_DIAMETER, 1, ...
    SOURCE_FREQ, SOUND_SPEED, DENSITY, x_vec, y_vec);

plot(x_vec, p_axial./max(p_axial))

xlabel('Axial distance [m]')
ylabel('Normalised pressure')
legend('AS k-Wave', 'O Neil soln')


