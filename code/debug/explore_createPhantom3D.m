% CREATEPHANTOM3D  Build full 3-D NIfTI phantoms for PRESTUS benchmarks.
%
% Writes final_tissues.nii.gz for benchmark subjects 001 (water/tissue) and
% 002 (water/skull/brain) into dedicated 3-D SimNIBS m2m folders.  The
% geometry is the rotationally symmetric extension of createPhantom.m:
% the same spherical-cap skull is rotated around the beam axis to produce a
% full bowl, and the same flat tissue interface is extruded in y.
%
% ── Array layout ────────────────────────────────────────────────────────────
%   dim 1 (x): lateral,  70 mm → 140 vox at 0.5 mm/vox
%   dim 2 (y): lateral,  70 mm → 140 vox at 0.5 mm/vox
%   dim 3 (z): axial,   120 mm → 240 vox at 0.5 mm/vox  (z=1: transducer)
%
% Coordinate grid: Xg/Yg/Zg = (1:N)*res mm, so voxel index 1 → 0.5 mm,
% voxel index 240 → 120 mm.
%
% ── NIfTI coordinate convention ─────────────────────────────────────────────
% write_phantom sets a diagonal sform (origin [0 0 0]) so that
%   world(i,j,k) = (i·res, j·res, k·res) mm.
% This matches the convention used by nifti_to_t1w for non-layered PRESTUS
% outputs, so the phantom overlays correctly in ITK-SNAP.
%
% ── Tissue labels (PRESTUS defaults) ────────────────────────────────────────
%   0 = water / coupling medium
%   1 = brain / soft tissue
%   4 = skull (cortical bone)
%
% ── Benchmark geometry ──────────────────────────────────────────────────────
%   Benchmark 1 (sub-101): flat water–tissue interface at 30 mm depth.
%     Focus at 64 mm depth lies 34 mm inside the soft-tissue slab.
%
%   Benchmark 2 (sub-102): spherical skull cap + brain.
%     Skull outer radius of curvature : 75.0 mm
%     Skull inner radius of curvature : 68.5 mm  (apex thickness = 6.5 mm)
%     Skull outer apex depth          : 30 mm
%     Centre of curvature             : 30 + 75 = 105 mm (on the brain side)
%     Brain fills everything inside the inner skull surface.
%
% See also: CREATEPHANTOM, PRESTUS_PARAMVAR_01_BENCHMARKS_SIM

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'))
rootpath = pwd;

pn.benchmark1 = fullfile(rootpath, '..', 'data', 'simnibs', 'm2m_sub-101');
pn.benchmark2 = fullfile(rootpath, '..', 'data', 'simnibs', 'm2m_sub-102');

for pf = {pn.benchmark1, pn.benchmark2}
    if ~exist(pf{1}, 'dir'); mkdir(pf{1}); end
end

%% grid

res = 0.5;                       % voxel size [mm]
nx  = round(70  / res);          % 140 — lateral x
ny  = round(70  / res);          % 140 — lateral y
nz  = round(120 / res);          % 240 — axial z (z = 1 is the transducer face)

% Physical coordinate (mm) of voxel index k along each axis: k * res.
% Voxel 1 → 0.5 mm, voxel 240 → 120 mm.
[Xg, Yg, Zg] = ndgrid((1:nx)*res, (1:ny)*res, (1:nz)*res);

%% Benchmark 1 — flat soft-tissue slab

data1 = zeros(nx, ny, nz, 'uint8');

% Tissue starts at 30 mm from the transducer face.
% Use strict > so the first tissue voxel is at Zg = 30.5 mm (voxel index 61),
% matching the 2-D phantom (z_tissue = round(30/res)+1 = 61).
z_interface_mm = 30;
data1(Zg > z_interface_mm) = 1;   % label 1 from 30.5 mm to domain end

write_phantom(data1, res, fullfile(pn.benchmark1, 'final_tissues.nii.gz'));
fprintf('Benchmark 1 3D (sub-101) written -> %s\n', pn.benchmark1);

%% Benchmark 2 — spherical skull cap + brain
%
% The skull is a spherical cap convex toward the transducer (anatomically
% realistic orientation).  The centre of curvature lies on the brain side:
%
%   cz_mm  =  apex_depth + r_outer  =  30 + 75  =  105 mm
%
% At the beam axis (Xg = cx_mm, Yg = cy_mm):
%   skull outer surface : Zg = 105 − 75  =  30.0 mm
%   skull inner surface : Zg = 105 − 68.5 =  36.5 mm
%   skull thickness     : 36.5 − 30.0     =   6.5 mm
%   brain onset         : Zg > 36.5 mm
%
% The posterior shell (Zg = 105 + 68.5 = 173.5 mm and beyond) falls outside
% the domain (nz = 120 mm), so no skull appears at the back of the volume.

cx_mm          = 35;      % lateral beam centre [mm]
cy_mm          = 35;
r_inner_mm     = 68.5;    % inner skull radius [mm]
r_outer_mm     = 75.0;    % outer skull radius [mm]
skull_apex_z_mm = 30;
cz_mm          = skull_apex_z_mm + r_outer_mm;   % 105 mm (on brain side)

% Distance from the centre of curvature at every grid point [mm].
dist_mm = sqrt((Xg - cx_mm).^2 + (Yg - cy_mm).^2 + (Zg - cz_mm).^2);

skull_mask = (dist_mm >= r_inner_mm) & (dist_mm <= r_outer_mm);  % shell
brain_mask = (dist_mm < r_inner_mm);                              % interior

data2 = zeros(nx, ny, nz, 'uint8');
data2(skull_mask) = 4;    % cortical bone
data2(brain_mask) = 1;    % soft tissue / brain

write_phantom(data2, res, fullfile(pn.benchmark2, 'final_tissues.nii.gz'));
fprintf('Benchmark 2 3D (sub-102) written -> %s\n', pn.benchmark2);

%% helper

function write_phantom(data, res, out_path)
% WRITE_PHANTOM  Save uint8 tissue array as a compressed NIfTI with sform.
%
% Sets a diagonal sform with origin [0, 0, 0] mm so that the phantom shares
% the same world-space convention as PRESTUS non-layered simulation outputs
% written by nifti_to_t1w:
%   world coordinate of voxel index (i, j, k):
%     x = i * res,   y = j * res,   z = k * res   [mm]
%
% This ensures that the tissue phantom overlays correctly on acoustic and
% thermal output maps when opened together in ITK-SNAP.

    % Diagonal affine: scale each axis by res, no translation.
    % affine3d expects the transpose of the standard 4×4 matrix
    % (MATLAB row-vector convention: [xin yin zin 1] * T = [xout yout zout 1]).
    T = diag([res res res 1]);

    tmp = strrep(out_path, '.nii.gz', '.nii');
    niftiwrite(double(data), tmp);            % write to get a valid header
    info                 = niftiinfo(tmp);
    info.PixelDimensions = repmat(res, 1, numel(info.PixelDimensions));
    info.Transform       = affine3d(T');
    niftiwrite(double(data), tmp, info);      % rewrite with corrected header
    gzip(tmp);
    delete(tmp);
end
