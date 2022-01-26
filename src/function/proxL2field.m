function out = proxL2field(y, gamma)
% Apply the proximation of the 2-norm to a vector field element by element
% (with respect to the 4th axis).
% proximity operator (proximity-operator.net)
    yNorm = sqrt(sum(y.^2, 4));
    s = max(0, 1 - gamma ./ yNorm);
    out = bsxfun(@times, y, s);
end
