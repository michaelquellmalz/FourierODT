classdef BackwardCGNE
    % Parameters for the backward model based on CGNE.
    properties
        CGNEopt = CGNEoptions
        enforceReal = true
        weightType = 'nonuniform'
        weights = []
    end
end

