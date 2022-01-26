function out = gradient(f)
% Compute the gradient of a matrix. 
% The partial derivaties are stored in the 4th axis.
    if ismatrix(f)
        out = cat(4, forwardDiffx(f), forwardDiffy(f));
    else
        out = cat(4, forwardDiffx(f), forwardDiffy(f), forwardDiffz(f));
    end
end

% forward finite differences (proximity-operator.net)
function out = forwardDiffx(f)
    out = [f(:,2:end,:) - f(:,1:end-1,:), zeros(size(f, 1), 1, size(f, 3))];
end

function out = forwardDiffy(f)
    out = [f(2:end,:,:) - f(1:end-1,:,:); zeros(1, size(f, 2), size(f, 3))];
end

function out = forwardDiffz(f)
    out = cat(3, f(:,:,2:end) - f(:,:,1:end-1), zeros(size(f, 1), size(f, 2), 1));
end