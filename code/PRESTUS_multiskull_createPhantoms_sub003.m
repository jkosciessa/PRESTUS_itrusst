% Create a 2D nifti image according to the desired tissue specification

% Here, we create skull media with cortical and trabecular subdivisions.
% This is naturally only sensible if a skull layer is present.
% We subdivide a 7mm thickness into 2mm cortical | 3mm trabecular | 2mm cortical
% For PRESTUS, cortical skull bone is index 7, trabecular is index 8

%% path management

% the following detects the path of this script (when run), which we can
% use to set up relative paths. in this example, PRESTUS is located in a
% directory called 'project/tools', and we would like to depoit the benchmark
% images into parallel folders project/data/simnibs/m2m-sub<xxx>

currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

%% Benchmarks 3xx

% This benchmark models a water / skull / soft tissue interface.
% It includes multiple thickness settings for the skull layer.
% The outer location of the skull is held constant.

pml = 20; % include a pml layer (on the axial dimension only)

thickness = 6.5;
i_thickness = 1;

distance_outer = 75;

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
distance_inner = distance_outer-thickness_current;

% Step 1: Initialize the matrix with global value 0 (dimensions: 120 rows (y), 70 columns (x))
data = zeros(120+pml, 70);

% Step 2: Generate coordinate grids for each pixel (assuming 1 mm per pixel)
[X, Y] = meshgrid(1:70, 1:120+pml);

% Step 3: Compute the center of curvature
center_x = 35;                  % x-coordinate (column)
center_y = 90 - 75;      % y-coordinate (row) => pml + 15

% Step 4: Calculate the Euclidean distance from the center
distance = sqrt((X - center_x).^2 + (Y - center_y).^2);

% Step 5: Assign values based on distance from center and y position
% Only fill below the center (Y >= center_y)

% ---- Skull subdivision: 2 mm cortical | 3 mm trabecular | 2 mm cortical ----
t_cort = 2;
t_trab = 2.5;

r_outer = distance_outer;
r_outer_cort_inner = r_outer - t_cort;                    % outer cortical inner radius
r_trab_inner       = r_outer_cort_inner - t_trab;         % trabecular inner radius
r_inner            = r_trab_inner - t_cort;               % should equal distance_inner if thickness=7

% Masks (concentric, no overlaps)
mask_outer_cort = (distance > r_outer_cort_inner) & (distance <= r_outer) & (Y >= center_y);
mask_trab      = (distance > r_trab_inner)       & (distance <= r_outer_cort_inner) & (Y >= center_y);
mask_inner_cort= (distance > r_inner)            & (distance <= r_trab_inner) & (Y >= center_y);

inner_mask     = (distance <= r_inner); % inner region (soft tissue)

% Labels (PRESTUS): cortical = 7, trabecular = 8
data(mask_outer_cort)  = 7;   % cortical (outer table)
data(mask_trab)        = 8;   % trabecular (diploë)
data(mask_inner_cort)  = 7;   % cortical (inner table)
data(inner_mask)       = 1;   % inner region (soft tissue)

% Step 6: Reshape data to 3D for NIfTI compatibility (singleton z-dimension)
data_3d = reshape(data, [pml + 120, 70, 1]);

% Invert image to match what we would like
data_3d = data_3d(end:-1:1,:,:);

% Step 7: Write the data to a NIfTI file with 1x1x1 mm voxel size
niftiwrite(data_3d, fullfile(pn.benchmark3, 'benchmark3.nii'));
gzip(fullfile(pn.benchmark3, 'benchmark3.nii'))
delete(fullfile(pn.benchmark3, 'benchmark3.nii'))
movefile(fullfile(pn.benchmark3, 'benchmark3.nii.gz'), ...
    fullfile(pn.benchmark3, 'final_tissues.nii.gz')) % rename to SimNIBS convention
