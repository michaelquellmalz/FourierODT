classdef BacktrackGLY < handle
    % Backtracking method of Goldstein, Li, Yuan (2015).
    properties (SetAccess = public)
        alpha0 = 0.95
        eta = 0.95
        c = 0.9
    end
    properties (SetAccess = protected)
        alpha
        tau
        sigma
        theta
        L
        G
        p
        d
        xNew
        xOld
        yNew
        yOld
        Dx
        Dy
        gradDx
    end
    methods
        function initParam(obj, tau, sigma, theta)
            obj.alpha = obj.alpha0;
            obj.tau = tau;
            obj.sigma = sigma;
            obj.theta = theta;
        end
        
        function initMap(obj, G, L)
            obj.G = G;
            obj.L = L;
        end
        
        function initPoint(obj, x, y)
            obj.xNew = x;
            obj.yNew = y;
        end
        
        function [tau, sigma] =update(obj, xNew, yNew)
            obj.updatePoints(xNew, yNew);
            obj.computeResidua;
            obj.backtracking;
            obj.balancing;
            tau = obj.tau;
            sigma = obj.sigma;
        end
    end
    
    methods (Access = protected)
        function updatePoints(obj, xNew, yNew)
            obj.xOld = obj.xNew;
            obj.xNew = xNew;
            obj.yOld = obj.yNew;
            obj.yNew = yNew;
        end
        
        function computeResidua(obj)
            obj.Dx = obj.xOld - obj.xNew;
            obj.Dy = obj.yOld - obj.yNew;
            obj.gradDx = obj.G.grad(obj.xOld) - obj.G.grad(obj.xNew);
            obj.p = obj.Dx / obj.tau - obj.L.adj(obj.Dy) - obj.gradDx; 
            obj.d = obj.Dy / obj.sigma - obj.theta * obj.L.op(obj.Dx);
        end
                
        function backtracking(obj)
            if obj.BTcond <= 0
                obj.tau = obj.tau / 2;
                obj.sigma = obj.sigma / 2;
            end
        end
        
        function B = BTcond(obj)
            B = obj.c / 2 / obj.tau * norm(obj.Dx(:))^2 ...
                + obj.c / 2 / obj.sigma * norm(obj.Dy(:))^2 ...
                - 2 * real(obj.Dy(:)' * reshape(obj.L.op(obj.Dx), [], 1)) ...
                - 2 * real(obj.Dx(:)' * obj.gradDx(:));
        end
        
        function balancing(obj)
            pN = norm(obj.p(:));
            dN = norm(obj.d(:));
            if pN >= 2 * dN
                obj.tau = obj.tau / (1 - obj.alpha);
                obj.sigma = obj.sigma * (1 - obj.alpha);
                obj.alpha = obj.alpha * obj.eta;
            elseif dN >= 2 * pN
                obj.tau = obj.tau * (1 - obj.alpha);
                obj.sigma = obj.sigma / (1 - obj.alpha);
                obj.alpha = obj.alpha * obj.eta;
            end
        end
    end
end

