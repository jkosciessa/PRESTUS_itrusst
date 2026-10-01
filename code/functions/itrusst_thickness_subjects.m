function [subject_ids, thickness_mm] = itrusst_thickness_subjects(thickness_mm)
% Map skull thickness (mm) to the subject IDs of the thickness phantoms.
%
%   integer mm: 300 + t          sub-301 = 1 mm   ... sub-312 = 12 mm
%   x.5 mm:     350 + floor(t)   sub-351 = 1.5 mm ... sub-361 = 11.5 mm
%
% Integer IDs are kept from the original 1 mm sweep, so existing phantoms
% and results (e.g. sub-303 = 3 mm) stay valid.
%
% Input:
%   thickness_mm - (optional) thickness values in mm, multiples of 0.5;
%                  default 1:0.5:12 (23 phantoms)
% Outputs:
%   subject_ids  - numeric subject IDs, same size as thickness_mm
%   thickness_mm - thickness values in mm (echoed back)

if nargin < 1
    thickness_mm = 1:0.5:12;
end

if any(mod(thickness_mm, 0.5) ~= 0) || any(thickness_mm < 1 | thickness_mm >= 50)
    error('Thickness must be a multiple of 0.5 mm in [1, 50).');
end

is_half     = mod(thickness_mm, 1) ~= 0;
subject_ids = 300 + floor(thickness_mm);
subject_ids(is_half) = 350 + floor(thickness_mm(is_half));
