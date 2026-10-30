function medium = itrusst_resolve_absorption(medium)
% Absorption coefficient of each tissue from its CSV specification.
%
% data/tissue_parameters.csv gives brain absorption as a fraction of attenuation
% (absorption_fraction; absorption_coeff empty) and skull and water absorption
% as a coefficient (absorption_coeff; absorption_fraction empty). Heat
% deposition uses the coefficient (thermal.absorption_model: coefficient), so
% a fraction is converted with the attenuation it is combined with:
%   absorption_coeff = absorption_fraction x alpha_coeff
% A brain attenuation sweep at a fixed fraction therefore scales absorption.
% A given coefficient is capped at the attenuation (absorption <= 100 %).
% Absorption and attenuation share the power law (alpha_power), so the
% coefficients compare directly. absorption_fraction is removed afterwards;
% PRESTUS does not read it with the coefficient model.
%
% Call after every override of alpha_coeff, absorption_coeff or
% absorption_fraction. Set both absorption fields of a tissue from the CSV:
% a finite absorption_coeff takes precedence, and the config defaults carry
% absorption_fraction values.
%
% Input/output:
%   medium - struct medium.(tissue).(property), as in parameters.medium_properties

tissues = fieldnames(medium);
for t = 1:numel(tissues)
    m = medium.(tissues{t});
    if ~isfield(m, 'alpha_coeff')
        continue
    end
    has_coeff    = isfield(m, 'absorption_coeff')    && ~isempty(m.absorption_coeff)    && ~isnan(m.absorption_coeff);
    has_fraction = isfield(m, 'absorption_fraction') && ~isempty(m.absorption_fraction) && ~isnan(m.absorption_fraction);
    if has_coeff
        m.absorption_coeff = min(m.absorption_coeff, m.alpha_coeff);
    elseif has_fraction
        if m.absorption_fraction < 0 || m.absorption_fraction > 1
            error('%s absorption_fraction = %g is outside [0, 1].', tissues{t}, m.absorption_fraction);
        end
        m.absorption_coeff = m.absorption_fraction * m.alpha_coeff;
    else
        error('%s has neither absorption_coeff nor absorption_fraction.', tissues{t});
    end
    if isfield(m, 'absorption_fraction')
        m = rmfield(m, 'absorption_fraction');
    end
    medium.(tissues{t}) = m;
end
