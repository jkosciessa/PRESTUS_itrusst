% CREATEPHANTOM  Build 2-D axisymmetric NIfTI phantoms for PRESTUS benchmarks.
%
% Writes final_tissues.nii.gz for benchmark subjects 001 (water/tissue) and
% 002 (water/skull/brain) into the SimNIBS m2m folders expected by PRESTUS.
%
% ── Array layout ────────────────────────────────────────────────────────────
%   dim 1 (x): lateral,  70 mm → 140 vox at 0.5 mm/vox  (full diameter)
%   dim 2    : singleton (size 1)
%   dim 3 (z): axial,   140 mm → 280 vox at 0.5 mm/vox  (z = 1: transducer)
%
% PRESTUS loads the file with squeeze(), yielding a [140 × 280] 2-D array.
% grid_axisymmetry then transposes and halves it to [Nz=280, Nr=70] before
% passing it to kspaceFirstOrderAS.
%
% ── NIfTI coordinate convention ─────────────────────────────────────────────
% write_phantom sets a diagonal sform so the singleton y-slice is centred at
% world y = (nx/2)*res mm, matching the y-centre of the 3-D output produced
% by convert_axisymmetric_to_3d + nifti_to_t1w.  This lets ITK-SNAP overlay
% the input tissue map on the simulated acoustic/thermal output.
%
% ── Tissue labels (PRESTUS defaults) ────────────────────────────────────────
%   0 = water / coupling medium
%   1 = brain / soft tissue
%   4 = skull (cortical bone)
%
% ── Benchmark geometry ──────────────────────────────────────────────────────
%   Benchmark 1 (sub-001): flat water–tissue interface at 30 mm depth.
%     Focus at 64 mm depth lies 34 mm inside the soft-tissue slab.
%
%   Benchmark 2 (sub-002): spherical skull cap + brain.
%     Skull outer radius of curvature : 75.0 mm
%     Skull inner radius of curvature : 68.5 mm  (apex thickness = 6.5 mm)
%     Skull outer apex depth          : 30 mm
%     Centre of curvature             : 30 + 75 = 105 mm (on the brain side)
%     Brain fills everything inside the inner skull surface.
%
% See also: CREATEPHANTOM3D, PRESTUS_PARAMVAR_01_BENCHMARKS_SIM

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr))
rootpath = pwd;

pn.benchmark1 = fullfile(rootpath, '..', 'data', 'simnibs', 'm2m_sub-001');
pn.benchmark2 = fullfile(rootpath, '..', 'data', 'simnibs', 'm2m_sub-002');

for pf = {pn.benchmark1, pn.benchmark2}
    if ~exist(pf{1}, 'dir'); mkdir(pf{1}); end
end

%% grid

res = 0.5;                     % voxel size [mm]
nx  = round(70  / res);        % 140 vox — full lateral diameter
nz  = round(140 / res);        % 280 vox — axial (z = 1 is the transducer face)
% PML is added outside the grid by k-Wave (PMLInside = false);
% these dimensions cover only the active simulation domain.

%% Benchmark 1 — water / soft-tissue interface

data1    = zeros(nx, 1, nz, 'uint8');
z_tissue = round(30 / res) + 1;   % first tissue voxel: 61 (= 30 mm depth)
data1(:, :, z_tissue:end) = 1;    % label 1 from 30 mm to domain end

write_phantom(data1, res, fullfile(pn.benchmark1, 'final_tissues.nii.gz'));
fprintf('Benchmark 1 written -> %s\n', pn.benchmark1);

%% Benchmark 2 — water / curved skull cap / brain
%
% The skull is a spherical cap convex toward the transducer (anatomically
% realistic orientation).  The centre of curvature lies on the brain side:
%
%   cz  =  apex_depth + r_outer  =  30/res + 150  =  60 + 150  =  210 vox
%        (= 105 mm, safely inside the 140-mm domain)
%
% At the beam axis (Xg = cx):
%   skull outer surface : Zg = cz − r_outer = 210 − 150 =  60 vox ( 30.0 mm)
%   skull inner surface : Zg = cz − r_inner = 210 − 137 =  73 vox ( 36.5 mm)
%   skull thickness     : 73 − 60 = 13 vox = 6.5 mm
%   brain onset         : Zg = 74 vox (37.0 mm)
%
% The skull shell does not wrap to the posterior side of the domain because
% cz + r_outer = 210 + 150 = 360 vox > nz = 280.

cx      = round(35 / res);     % lateral beam axis: vox 70 (= 35 mm)
r_inner = 68.5 / res;          % inner skull radius: 137 vox
r_outer = 75.0 / res;          % outer skull radius: 150 vox
cz      = 30 / res + r_outer;  % centre of curvature: 210 vox (= 105 mm)

% Distance from the centre of curvature at every grid point.
% ndgrid returns 2-D arrays [nx × nz]; dim 2 (singleton) is handled by
% indexing the 3-D data array separately below.
[Xg, Zg] = ndgrid(1:nx, 1:nz);
dist      = sqrt((Xg - cx).^2 + (Zg - cz).^2);

skull_mask = (dist >= r_inner) & (dist <= r_outer);  % spherical shell
brain_mask = (dist < r_inner);                        % interior of shell

data2 = zeros(nx, 1, nz, 'uint8');
data2(skull_mask) = 4;    % cortical bone
data2(brain_mask) = 1;    % soft tissue / brain

write_phantom(data2, res, fullfile(pn.benchmark2, 'final_tissues.nii.gz'));
fprintf('Benchmark 2 written -> %s\n', pn.benchmark2);

%% helper

function write_phantom(data, res, out_path)
% WRITE_PHANTOM  Save uint8 tissue array as a compressed NIfTI with sform.
%
% The NIfTI sform is set so that the singleton y-slice (dim 2 = 1 voxel)
% aligns spatially with the centre cross-section of the 3-D output produced
% by PRESTUS after radial expansion (convert_axisymmetric_to_3d).
%
% Convention (matches nifti_to_t1w non-layered path):
%   world coordinate of voxel index (i, j, k):
%     x = i * res + ox,   ox = 0
%     y = j * res + oy,   oy = (nx/2 − 1) * res   ← centres the slice
%     z = k * res + oz,   oz = 0
%
% The 3-D PRESTUS output has nx voxels in y with origin oy = 0, so its
% y-centre is at (nx/2) * res.  Setting oy = (nx/2 − 1)*res means the
% phantom's single slice (j = 1) maps to world y = 1*res + oy = (nx/2)*res,
% matching the output centre slice exactly.

    nx  = size(data, 1);
    oy  = (nx/2 - 1) * res;   % y-offset [mm] — aligns with 3-D output centre

    % Diagonal affine: scale by res along each axis, shift y by oy.
    % affine3d expects the transpose of the standard 4×4 matrix.
    T        = diag([res res res 1]);
    T(2, 4)  = oy;

    tmp = strrep(out_path, '.nii.gz', '.nii');
    niftiwrite(double(data), tmp);            % write to get a valid header
    info                 = niftiinfo(tmp);
    info.PixelDimensions = repmat(res, 1, numel(info.PixelDimensions));
    info.Transform       = affine3d(T');
    niftiwrite(double(data), tmp, info);      % rewrite with corrected header
    gzip(tmp);
    delete(tmp);
end
