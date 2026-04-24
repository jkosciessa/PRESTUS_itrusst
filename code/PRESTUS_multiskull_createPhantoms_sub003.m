% Create a 2D nifti image according to the desired tissue specification
% Resolution: 0.5 mm/voxel

% Here, we create skull media with cortical and trabecular subdivisions.
% This is naturally only sensible if a skull layer is present.
% We subdivide a 6.5mm thickness into 2mm cortical | 2.5mm trabecular | 2mm cortical
% For PRESTUS, cortical skull bone is index 7, trabecular is index 8

%% path management

currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

%% Benchmarks 3xx

% This benchmark models a water / skull / soft tissue interface.
% It includes multiple thickness settings for the skull layer.
% The outer location of the skull is held constant.

res = 0.5; % voxel size [mm]
% PML is added outside the grid by k-Wave (PMLInside=false).
% Phantom contains only the active domain — no embedded PML voxels.

thickness = 6.5;
i_thickness = 1;

distance_outer = 75; % [mm]

% Provide a subject label: sub-003
pn.benchmark3 = fullfile(rootpath, 'data', 'simnibs', 'm2m_sub-003');
if isempty(pn.benchmark3)
    pn.benchmark3 = rootpath;
else
    if ~exist(pn.benchmark3)
        mkdir(pn.benchmark3);
    end
end

thickness_current = thickness(i_thickness);
distance_inner = distance_outer - thickness_current;

% Active domain size (= full grid, PML is external)
nx = round(70  / res);   % 140 voxels = 70 mm lateral
ny = round(120 / res);   % 240 voxels = 120 mm axial

% Step 1: Initialize the matrix with global value 0
data = zeros(ny, nx);

% Step 2: Generate coordinate grids
[X, Y] = meshgrid(1:nx, 1:ny);

% Step 3: Compute the center of curvature within the active domain (no PML offset)
center_x = round(35 / res);
center_y = round(15 / res);

% Step 4: Calculate the Euclidean distance from the center
distance = sqrt((X - center_x).^2 + (Y - center_y).^2);

% Step 5: Assign values based on distance from center and y position

% ---- Skull subdivision: 2 mm cortical | 2.5 mm trabecular | 2 mm cortical ----
t_cort = 2;
t_trab = 2.5;

r_outer            = distance_outer / res;
r_outer_cort_inner = r_outer - t_cort / res;             % outer cortical inner radius
r_trab_inner       = r_outer_cort_inner - t_trab / res;  % trabecular inner radius
r_inner            = r_trab_inner - t_cort / res;        % should equal distance_inner/res if thickness=6.5

% Masks (concentric, no overlaps)
mask_outer_cort = (distance > r_outer_cort_inner) & (distance <= r_outer)            & (Y >= center_y);
mask_trab       = (distance > r_trab_inner)       & (distance <= r_outer_cort_inner) & (Y >= center_y);
mask_inner_cort = (distance > r_inner)            & (distance <= r_trab_inner)       & (Y >= center_y);
inner_mask      = (distance <= r_inner);

% Labels (PRESTUS): cortical = 7, trabecular = 8
data(mask_outer_cort) = 7;   % cortical (outer table)
data(mask_trab)       = 8;   % trabecular (diploë)
data(mask_inner_cort) = 7;   % cortical (inner table)
data(inner_mask)      = 1;   % inner region (soft tissue)

% Step 6: Reshape data to 3D for NIfTI compatibility: [nx x 1 x ny]
% Singleton in the second dimension so external tools can overlay this
% 2D phantom slice with 3D simulation results [nx x nx x ny].
data_3d = reshape(data.', [nx, 1, ny]);

% Invert image to match what we would like
data_3d = data_3d(:,:,end:-1:1);

% Step 7: Write the data to a NIfTI file
% Note: physical resolution is set via grid.resolution_mm in the PRESTUS config
niftiwrite(data_3d, fullfile(pn.benchmark3, 'benchmark3.nii'));
gzip(fullfile(pn.benchmark3, 'benchmark3.nii'))
delete(fullfile(pn.benchmark3, 'benchmark3.nii'))
movefile(fullfile(pn.benchmark3, 'benchmark3.nii.gz'), ...
    fullfile(pn.benchmark3, 'final_tissues.nii.gz')) % rename to SimNIBS convention
