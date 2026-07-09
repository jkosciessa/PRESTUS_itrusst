% CREATEPHANTOM  Build 2-D axisymmetric NIfTI phantoms for PRESTUS benchmarks.
%
% Writes final_tissues.nii.gz for benchmark subjects 001 (water/tissue) and
% 002 (water/skull/brain) into the SimNIBS m2m folders expected by PRESTUS.
%
% ── Array layout ────────────────────────────────────────────────────────────
%   dim 1 (x): lateral,  70 mm → 140 vox at 0.5 mm/vox  (full diameter)
%   dim 2 (z): axial,   120 mm → 240 vox at 0.5 mm/vox  (z = 1: transducer)
%
% PRESTUS loads the file with squeeze(), yielding a [140 × 240] 2-D array.
% grid_axisymmetry then transposes and halves it to [Nz=240, Nr=70] before
% passing it to kspaceFirstOrderAS.
%
% ── NIfTI coordinate convention ─────────────────────────────────────────────
% Phantoms are stored as true 2-D arrays (no singleton dimension).
% nifti_to_t1w copies phantom_header.Transform to all simulation outputs, so
% phantom and outputs share the same affine and align directly in ITK-SNAP.
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
nz  = round(120 / res);        % 240 vox — axial (z = 1 is the transducer face)
% PML is added outside the grid by k-Wave (PMLInside = false);
% these dimensions cover only the active simulation domain.

%% Benchmark 1 — water / soft-tissue interface

data1    = zeros(nx, nz, 'uint8');
z_tissue = round(30 / res) + 1;   % first tissue voxel: 61 (= 30 mm depth)
data1(:, z_tissue:end) = 1;       % label 1 from 30 mm to domain end

write_phantom(data1, res, fullfile(pn.benchmark1, 'final_tissues.nii.gz'));
fprintf('Benchmark 1 written -> %s\n', pn.benchmark1);

%% Benchmark 2 — water / curved skull cap / brain
%
% The skull is a spherical cap convex toward the transducer (anatomically
% realistic orientation).  The centre of curvature lies on the brain side:
%
%   cz  =  apex_depth + r_outer  =  30/res + 150  =  60 + 150  =  210 vox
%        (= 105 mm, safely inside the 120-mm domain)
%
% At the beam axis (Xg = cx):
%   skull outer surface : Zg = cz − r_outer = 210 − 150 =  60 vox ( 30.0 mm)
%   skull inner surface : Zg = cz − r_inner = 210 − 137 =  73 vox ( 36.5 mm)
%   skull thickness     : 73 − 60 = 13 vox = 6.5 mm
%   brain onset         : Zg = 74 vox (37.0 mm)
%
% The skull shell does not wrap to the posterior side of the domain because
% cz + r_outer = 210 + 150 = 360 vox > nz = 240.

cx      = round(35 / res);     % lateral beam axis: vox 70 (= 35 mm)
r_inner = 68.5 / res;          % inner skull radius: 137 vox
r_outer = 75.0 / res;          % outer skull radius: 150 vox
cz      = 30 / res + r_outer;  % centre of curvature: 210 vox (= 105 mm)

[Xg, Zg] = ndgrid(1:nx, 1:nz);
dist      = sqrt((Xg - cx).^2 + (Zg - cz).^2);

skull_mask = (dist >= r_inner) & (dist <= r_outer);  % spherical shell
brain_mask = (dist < r_inner);                        % interior of shell

data2 = zeros(nx, nz, 'uint8');
data2(skull_mask) = 4;    % cortical bone
data2(brain_mask) = 1;    % soft tissue / brain

write_phantom(data2, res, fullfile(pn.benchmark2, 'final_tissues.nii.gz'));
fprintf('Benchmark 2 written -> %s\n', pn.benchmark2);

%% helper

function write_phantom(data, res, out_path)
% WRITE_PHANTOM  Save uint8 tissue array as a compressed NIfTI.
%
% Writes a diagonal-sform NIfTI with isotropic voxel size res. No singleton
% dimension: data is [nx × nz] where dim 1 = lateral (x) and dim 2 = axial.
% nifti_to_t1w copies this header to all simulation outputs, so the phantom
% and outputs share the same affine and align directly in ITK-SNAP.

    tmp = strrep(out_path, '.nii.gz', '.nii');
    niftiwrite(double(data), tmp);
    info                 = niftiinfo(tmp);
    info.PixelDimensions = repmat(res, 1, numel(info.PixelDimensions));
    niftiwrite(double(data), tmp, info);
    gzip(tmp);
    delete(tmp);
end
