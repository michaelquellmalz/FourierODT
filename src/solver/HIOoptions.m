classdef HIOoptions
    properties
        updateMethod = 'HIO'    % HIO - Hybrid-Input-Output
                                % BIO - Basic-Input-Output
                                % OO - Output-Output
        beta = 0.7
        maxIter = 100
        tol = 1e-8
        enforceReal = true
        enforceNonneg = true
        supportRadius = 6
        condTol = 1e-10         % Tolerence for real and support condition.
        f0 = []
    end
end

