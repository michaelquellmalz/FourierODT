classdef BackwardFBPD
    % Parameters for the backward model based on FBPD.
    properties
        FBPDopt = FBPDoptions
        enforceReal = true
        enforceNonneg = true
        lambda = 0.5
        weightType = 'uniform'
        weights = []
        resume = false
        state = []
    end
end

