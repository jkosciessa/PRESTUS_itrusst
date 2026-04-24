% Create a 2D nifti image according to the desired tissue specification

% Here, we create skull media with different thickness.
% This is naturally only sensible if a skull layer is present.

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

res = 0.5; % voxel size [mm]
% PML is added outside the grid by k-Wave (PMLInside=false).
% Phantoms contain only the active domain — no embedded PML voxels.

thickness = 1:1:12;
distance_outer = 75; % [mm]

% Active domain size (= full grid, PML is external)
nx = round(70  / res);   % 140 voxels = 70 mm lateral
ny = round(120 / res);   % 240 voxels = 120 mm axial

for i_thickness = 1:numel(thickness)

    % Provide each thickness a subject label sub-3xx
    pn.benchmark3 = fullfile(rootpath, 'data', 'simnibs', sprintf('m2m_sub-3%02d', i_thickness));
    if isempty(pn.benchmark3)
        pn.benchmark3 = rootpath;
    else
        if ~exist(pn.benchmark3)
            mkdir(pn.benchmark3);
        end
    end

    thickness_current = thickness(i_thickness);
    distance_inner = distance_outer - thickness_current;

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
    % Only fill below the center (Y >= center_y)
    annulus_mask = (distance > distance_inner / res) & (distance <= distance_outer / res) & (Y >= center_y);
    inner_mask   = (distance <= distance_inner / res);

    data(annulus_mask) = 4;    % Annular layer
    data(inner_mask)   = 1;    % Inner region

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
end