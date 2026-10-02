% Create a 2D nifti image according to the desired tissue specification
%
% Creates skull media with different thickness (1–12 mm sweep, 0.5 mm steps).
% Subject IDs follow itrusst_thickness_subjects.m: sub-301..312 = 1..12 mm,
% sub-351..361 = 1.5..11.5 mm.
%
% ── Array layout ─────────────────────────────────────────────────────────────
%   dim 1 (x): lateral,  70 mm → 140 vox at 0.5 mm/vox  (full diameter)
%   dim 2 (z): axial,   120 mm → 240 vox at 0.5 mm/vox  (z = 1: transducer)
%
% Mirrors the creation principle of createPhantom.m:
%   no axis flip, no singleton dimension, same write_phantom convention.
%
% ── Benchmark geometry ───────────────────────────────────────────────────────
%   Outer skull surface held constant at r_outer = 75 mm, apex at 30 mm depth.
%   Inner surface varies: distance_inner = 75 - thickness.
%   cz + r_outer = 210 + 150 = 360 vox > nz = 240 → sphere truncated by domain,
%   no half-cap guard needed.

%% path management

currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

%% grid

res = 0.5;               % voxel size [mm]
nx  = round(70  / res);  % 140 vox — full lateral diameter
nz  = round(120 / res);  % 240 vox — axial (z = 1: transducer face)

%% geometry constants

skull_apex_z_mm = 30;
distance_outer  = 75;   % outer radius of curvature [mm]

cx = round(35 / res);                        % lateral beam axis: vox 70
r_outer = distance_outer / res;              % 150 vox
cz = skull_apex_z_mm / res + r_outer;       % centre of curvature: vox 210 (= 105 mm)

[Xg, Zg] = ndgrid(1:nx, 1:nz);

%% thickness sweep

addpath(fullfile(rootpath, 'code', 'functions'));
[subject_ids, thickness] = itrusst_thickness_subjects(1:0.5:12);   % [mm]

for i_thickness = 1:numel(thickness)

    pn.benchmark3 = fullfile(rootpath, 'data', 'simnibs', sprintf('m2m_sub-%03d', subject_ids(i_thickness)));
    if ~exist(pn.benchmark3, 'dir')
        mkdir(pn.benchmark3);
    end

    thickness_current = thickness(i_thickness);
    distance_inner    = distance_outer - thickness_current;
    r_inner           = distance_inner / res;

    dist = sqrt((Xg - cx).^2 + (Zg - cz).^2);

    % No half-cap guard — domain truncates the sphere naturally
    annulus_mask = (dist > r_inner) & (dist <= r_outer);
    inner_mask   = (dist <= r_inner);

    data = zeros(nx, nz, 'uint8');
    data(annulus_mask) = 4;   % skull (cortical bone)
    data(inner_mask)   = 1;   % soft tissue / brain

    write_phantom(data, res, fullfile(pn.benchmark3, 'final_tissues.nii.gz'));
    fprintf('Benchmark 3 thickness %.1f mm (sub-%03d) written -> %s\n', ...
        thickness_current, subject_ids(i_thickness), pn.benchmark3);

end

%% helper

function write_phantom(data, res, out_path)
% WRITE_PHANTOM  Save uint8 tissue array as a compressed NIfTI.
% No singleton dimension: data is [nx × nz] (lateral × axial).
% Matches createPhantom.m convention.

    tmp = strrep(out_path, '.nii.gz', '.nii');
    niftiwrite(double(data), tmp);
    info                 = niftiinfo(tmp);
    info.PixelDimensions = repmat(res, 1, numel(info.PixelDimensions));
    niftiwrite(double(data), tmp, info);
    gzip(tmp);
    delete(tmp);
end
