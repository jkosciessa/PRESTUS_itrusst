% Create a 2D nifti image according to the desired tissue specification
% Resolution: 0.25 mm/voxel
%
% Creates skull media with cortical and trabecular subdivisions.
% Skull layer (6.5 mm total): 2 mm cortical | 2.5 mm trabecular | 2 mm cortical
% PRESTUS tissue labels: cortical = 7, trabecular = 8
%
% ── Array layout ─────────────────────────────────────────────────────────────
%   dim 1 (x): lateral,  70 mm → 280 vox at 0.25 mm/vox  (full diameter)
%   dim 2 (z): axial,   120 mm → 480 vox at 0.25 mm/vox  (z = 1: transducer)
%
% Mirrors the creation principle of createPhantom.m / createPhantom025mm.m:
%   no axis flip, no singleton dimension, same write_phantom convention.
%
% ── Benchmark geometry ───────────────────────────────────────────────────────
%   Skull outer radius of curvature : 75.0 mm
%   Skull inner radius of curvature : 68.5 mm  (apex thickness = 6.5 mm)
%   Skull outer apex depth          : 30 mm
%   Centre of curvature             : 30 + 75 = 105 mm (on the brain side)
%   cz + r_outer = 420 + 300 = 720 vox > nz = 480 → sphere truncated by domain,
%   no half-cap guard needed.

%% path management

currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

%% grid

res = 0.25;              % voxel size [mm]
nx  = round(70  / res);  % 280 vox — full lateral diameter
nz  = round(120 / res);  % 480 vox — axial (z = 1: transducer face)

%% subject folder

pn.benchmark3 = fullfile(rootpath, 'data', 'simnibs', 'm2m_sub-003');
if ~exist(pn.benchmark3, 'dir')
    mkdir(pn.benchmark3);
end

%% geometry

skull_apex_z_mm = 30;
thickness       = 6.5;   % total skull thickness [mm]
distance_outer  = 75;    % outer radius of curvature [mm]
distance_inner  = distance_outer - thickness;

cx      = round(35 / res);                       % lateral beam axis: vox 140
r_outer = distance_outer / res;                  % 300 vox
r_inner = distance_inner / res;                  % 274 vox
cz      = skull_apex_z_mm / res + r_outer;       % centre of curvature: vox 420 (= 105 mm)

% Skull subdivision
t_cort = 2;    % cortical layer thickness [mm]
t_trab = 2.5;  % trabecular layer thickness [mm]

r_outer_cort_inner = r_outer - t_cort / res;              % 292 vox
r_trab_inner       = r_outer_cort_inner - t_trab / res;   % 282 vox
% r_inner should equal r_trab_inner - t_cort/res = 274 (verified)

[Xg, Zg] = ndgrid(1:nx, 1:nz);
dist = sqrt((Xg - cx).^2 + (Zg - cz).^2);

% Concentric shell masks (no half-cap guard — domain truncates the sphere)
mask_outer_cort = (dist > r_outer_cort_inner) & (dist <= r_outer);
mask_trab       = (dist > r_trab_inner)       & (dist <= r_outer_cort_inner);
mask_inner_cort = (dist > r_inner)            & (dist <= r_trab_inner);
inner_mask      = (dist <= r_inner);

data = zeros(nx, nz, 'uint8');
data(mask_outer_cort) = 7;   % cortical (outer table)
data(mask_trab)       = 8;   % trabecular (diploë)
data(mask_inner_cort) = 7;   % cortical (inner table)
data(inner_mask)      = 1;   % soft tissue / brain

write_phantom(data, res, fullfile(pn.benchmark3, 'final_tissues.nii.gz'));
fprintf('Benchmark 3 multiskull (sub-003) written -> %s\n', pn.benchmark3);

%% helper

function write_phantom(data, res, out_path)
% WRITE_PHANTOM  Save uint8 tissue array as a compressed NIfTI.
% No singleton dimension: data is [nx × nz] (lateral × axial).
% Matches createPhantom.m / createPhantom025mm.m convention.

    tmp = strrep(out_path, '.nii.gz', '.nii');
    niftiwrite(double(data), tmp);
    info                 = niftiinfo(tmp);
    info.PixelDimensions = repmat(res, 1, numel(info.PixelDimensions));
    niftiwrite(double(data), tmp, info);
    gzip(tmp);
    delete(tmp);
end
