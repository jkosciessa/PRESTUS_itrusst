% Create a 2D nifti image according to the desired tissue specification
% Resolution: 0.5 mm/voxel

%% path management

% the following detects the path of this script (when run), which we can
% use to set up relative paths. in this example, PRESTUS is located in a
% directory called 'project/tools', and we would like to depoit the benchmark
% images into parallel folders project/data/simnibs/m2m-sub<xxx>

currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr))
rootpath = pwd;

pn.benchmark1 = fullfile(rootpath, '..', 'data', 'simnibs', 'm2m_sub-001'); % specify the path to a SimNIBS m2m folder (optional)
pn.benchmark2 = fullfile(rootpath, '..', 'data', 'simnibs', 'm2m_sub-002'); % specify second benchmark as a different subject (optional)

% if no path is provided, the folder containing PRESTUS is specified by default
% if paths are specified and do not yet exist, create the folders
if isempty(pn.benchmark1)
    pn.benchmark1 = rootpath;
else
    if ~exist(pn.benchmark1)
        mkdir(pn.benchmark1);
    end
end
if isempty(pn.benchmark2)
    pn.benchmark2 = rootpath;
else
    if ~exist(pn.benchmark2)
        mkdir(pn.benchmark2);
    end
end

res = 0.5; % voxel size [mm]
% PML is added outside the grid by k-Wave (PMLInside=false).
% Phantoms contain only the active domain — no embedded PML voxels.

%% Benchmark 1

% This benchmark models a water / soft tissue interface.

% Active domain size (= full grid, PML is external)
nx = round(70  / res);   % 140 voxels = 70 mm lateral
ny = round(120 / res);   % 240 voxels = 120 mm axial

% Initialize the matrix with global value 0
data = zeros(ny, nx);

% Generate coordinate grids
[X, Y] = meshgrid(1:nx, 1:ny);

% Fill soft tissue region (starts at 30 mm into active domain)
row_start = round(30 / res) + 1;
data(row_start : ny, :) = 1;

% Reshape data to 3D for NIfTI compatibility: [nx x 1 x ny]
% Singleton in the second dimension so external tools can overlay this
% 2D phantom slice with 3D simulation results [nx x nx x ny].
data_3d = reshape(data.', [nx, 1, ny]);

% Write the data to a NIfTI file
% Note: physical resolution is set via grid.resolution_mm in the PRESTUS config
niftiwrite(data_3d, fullfile(pn.benchmark1, 'benchmark1.nii'));
gzip(fullfile(pn.benchmark1, 'benchmark1.nii')) % Note: PRESTUS expects compressed images
delete(fullfile(pn.benchmark1, 'benchmark1.nii')) % remove uncompressed image
movefile(fullfile(pn.benchmark1, 'benchmark1.nii.gz'), ...
    fullfile(pn.benchmark1, 'final_tissues.nii.gz')) % rename to SimNIBS convention

% Note: the output file should be called final_tissues.nii.gz and be
% located in a "subject-"/benchmark specific simnibs m2m_sub-xxx folder

%% Benchmark 2

% This benchmark models a water / skull / soft tissue interface.

% Initialize the active grid
data = zeros(ny, nx);

% Generate coordinate grids
[X, Y] = meshgrid(1:nx, 1:ny);

% Center of curvature within the active domain (no PML offset)
center_x = round(35 / res);
center_y = round(15 / res);

% Euclidean distance from the center
distance = sqrt((X - center_x).^2 + (Y - center_y).^2);

% Skull annulus: inner radius 68.5 mm, outer radius 75 mm
r_inner = 68.5 / res;
r_outer = 75.0 / res;
annulus_mask = (distance > r_inner) & (distance <= r_outer) & (Y >= center_y);
inner_mask   = (distance <= r_inner);

data(annulus_mask) = 4;    % Annular layer
data(inner_mask)   = 1;    % Inner region

% Reshape data to 3D for NIfTI compatibility: [nx x 1 x ny]
data_3d = reshape(data.', [nx, 1, ny]);

% Invert image to match what we would like
data_3d = data_3d(:,:,end:-1:1);

% Write the data to a NIfTI file
% Note: physical resolution is set via grid.resolution_mm in the PRESTUS config
niftiwrite(data_3d, fullfile(pn.benchmark2, 'benchmark2.nii'));
gzip(fullfile(pn.benchmark2, 'benchmark2.nii'))
delete(fullfile(pn.benchmark2, 'benchmark2.nii'))
movefile(fullfile(pn.benchmark2, 'benchmark2.nii.gz'), ...
    fullfile(pn.benchmark2, 'final_tissues.nii.gz')) % rename to SimNIBS convention
