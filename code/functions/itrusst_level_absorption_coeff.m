function T = itrusst_level_absorption_coeff(T)
% Absorption coefficient per row of tissue_parameters.csv, for range figures.
%
% Rows with an absorption_fraction (brain) get absorption_coeff =
% absorption_fraction x alpha_coeff of the same level (bottom/mid/top); rows
% with a coefficient (skull, water) keep it. Simulations combine the fraction
% with the attenuation in use instead (itrusst_resolve_absorption.m).
%
% Input/output:
%   T - table read from data/tissue_parameters.csv

from_fraction = isnan(T.absorption_coeff) & ~isnan(T.absorption_fraction);
T.absorption_coeff(from_fraction) = T.absorption_fraction(from_fraction) .* T.alpha_coeff(from_fraction);
