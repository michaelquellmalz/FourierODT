classdef BacktrackYH < handle
    % Backtracking method of Yokota, Hontani (2015).
    properties (SetAccess = public)
        rho = 0.005
        c = 0.9
        beta = 1.5
        zeta = 0.25
    end
    properties (SetAccess = protected)
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
            %obj.balancing;
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
            wp = real(obj.Dx(:)' * obj.p(:)) / norm(obj.Dx(:)) / norm(obj.p(:));
            wd = real(obj.Dy(:)' * obj.d(:)) / norm(obj.Dy(:)) / norm(obj.d(:));
            if wp > obj.c
                obj.tau = obj.beta * obj.tau;
            elseif wp < 0
                obj.tau = obj.zeta * obj.tau;
            end
            if wd > obj.c
                obj.sigma = obj.beta * obj.sigma;
            elseif wd < 0
                obj.sigma = obj.zeta * obj.sigma;
            end
        end
        
        function balancing(obj)
            R = (norm(obj.p(:)) / norm(obj.d(:)))^obj.rho;
            obj.tau = obj.tau * R;
            obj.sigma = obj.sigma / R;
        end
    end
end

