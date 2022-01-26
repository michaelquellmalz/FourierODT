classdef TVoptions
    % Parameters for the backward model based on FBPD.
    properties
        FBPDopt = FBPDoptions
        enforceReal = true
        enforceNonneg = true
        lambda = 0.5
    end
end

