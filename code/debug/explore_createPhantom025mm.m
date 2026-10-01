% CREATEPHANTOM025MM  Build 0.25 mm 2-D axisymmetric NIfTI phantoms for PRESTUS benchmarks.
%
% Writes final_tissues.nii.gz for benchmark subjects 201 (water/tissue) and
% 202 (water/skull/brain) into the SimNIBS m2m folders expected by PRESTUS.
% This is identical in geometry to createPhantom.m but uses 0.25 mm/vox
% resolution (12 PPW at 500 kHz) to assess the effect of resolution on the
% jagged on-axis pressure profiles observed at 6 PPW (0.5 mm/vox).
%
% ── Array layout ────────────────────────────────────────────────────────────
%   dim 1 (x): lateral,  70 mm → 280 vox at 0.25 mm/vox  (full diameter)
%   dim 2 (z): axial,   140 mm → 560 vox at 0.25 mm/vox  (z = 1: transducer)
%
% ── Geometry ────────────────────────────────────────────────────────────────
% Identical to createPhantom.m:
%   Benchmark 1 (sub-201): flat water–tissue interface at 30 mm depth.
%   Benchmark 2 (sub-202): spherical skull cap (r_outer=75 mm, r_inner=68.5 mm,
%     apex at 30 mm depth) + brain.
%
% See also: CREATEPHANTOM, CREATEPHANTOM3D, PRESTUS_PARAMVAR_01_BENCHMARKS_025MM_SIM

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'))
rootpath = pwd;

pn.benchmark1 = fullfile(rootpath, '..', 'data', 'simnibs', 'm2m_sub-201');
pn.benchmark2 = fullfile(rootpath, '..', 'data', 'simnibs', 'm2m_sub-202');

for pf = {pn.benchmark1, pn.benchmark2}
    if ~exist(pf{1}, 'dir'); mkdir(pf{1}); end
end

%% grid

res = 0.25;                    % voxel size [mm]
nx  = round(70  / res);        % 280 vox — full lateral diameter
nz  = round(140 / res);        % 560 vox — axial (z = 1 is the transducer face)

%% Benchmark 1 — water / soft-tissue interface

data1    = zeros(nx, nz, 'uint8');
z_tissue = round(30 / res) + 1;   % first tissue voxel: 121 (= 30 mm depth)
data1(:, z_tissue:end) = 1;

write_phantom(data1, res, fullfile(pn.benchmark1, 'final_tissues.nii.gz'));
fprintf('Benchmark 1 0.25 mm (sub-201) written -> %s\n', pn.benchmark1);

%% Benchmark 2 — water / curved skull cap / brain

cx      = round(35 / res);     % lateral beam axis: vox 140 (= 35 mm)
r_inner = 68.5 / res;          % inner skull radius: 274 vox
r_outer = 75.0 / res;          % outer skull radius: 300 vox
cz      = 30 / res + r_outer;  % centre of curvature: 420 vox (= 105 mm)

[Xg, Zg] = ndgrid(1:nx, 1:nz);
dist      = sqrt((Xg - cx).^2 + (Zg - cz).^2);

skull_mask = (dist >= r_inner) & (dist <= r_outer);
brain_mask = (dist < r_inner);

data2 = zeros(nx, nz, 'uint8');
data2(skull_mask) = 4;
data2(brain_mask) = 1;

write_phantom(data2, res, fullfile(pn.benchmark2, 'final_tissues.nii.gz'));
fprintf('Benchmark 2 0.25 mm (sub-202) written -> %s\n', pn.benchmark2);

%% helper

function write_phantom(data, res, out_path)
% WRITE_PHANTOM  Save uint8 tissue array as a compressed NIfTI.
% No singleton dimension: data is [nx × nz] (lateral × axial).

    tmp = strrep(out_path, '.nii.gz', '.nii');
    niftiwrite(double(data), tmp);
    info                 = niftiinfo(tmp);
    info.PixelDimensions = repmat(res, 1, numel(info.PixelDimensions));
    niftiwrite(double(data), tmp, info);
    gzip(tmp);
    delete(tmp);
end
