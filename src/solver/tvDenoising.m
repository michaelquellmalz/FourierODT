function f = tvDenoising(f0, opt)
    F = buildF(opt);
    G = buildG(f0);
    H = buildH(opt);
    L = buildL;
    PD = FBPD(F, G, H, L, opt.FBPDopt);
    f = PD.minimize(f0);
end

function F = buildF(opt)
    if opt.enforceNonneg
        F.prox = @(f, tau) max(0, real(f));
    elseif opt.enforceReal
        F.prox = @(f, tau) real(f);
    else
        F.prox = @(f, tau) f;
    end
end

function G = buildG(f0)
    G.grad = @(f) f - f0;
end

function H = buildH(opt)
    lambda = opt.lambda;
    H.prox = @(y,gamma) proxL2field(y, gamma*lambda);
end

function L = buildL()
    L.op = @(f) gradient(f);
    L.adj = @(y) divergence(y);
end