function out = divergence(y)
% Compute divergence of a vector field (with respect to 4th axis).
    if size(y, 3) == 1
        out = backwardDiffx(y(:,:,:,1)) + backwardDiffy(y(:,:,:,2));
    else
        out = backwardDiffx(y(:,:,:,1)) + backwardDiffy(y(:,:,:,2)) ...
            + backwardDiffz(y(:,:,:,3));        
    end
end

% backward finite differences (proximity-operator.net)
function out = backwardDiffx(f) 
    out = [-f(:,1,:), f(:,1:end-2,:) - f(:,2:end-1,:), f(:,end-1,:)];
end

function out = backwardDiffy(f)
    out = [-f(1,:,:); f(1:end-2,:,:) - f(2:end-1,:,:); f(end-1,:,:)];
end

function out = backwardDiffz(f)
    out = cat(3, -f(:,:,1), f(:,:,1:end-2) - f(:,:,2:end-1), f(:,:,end-1));
end