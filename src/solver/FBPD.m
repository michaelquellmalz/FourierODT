classdef FBPD < handle
    properties (SetAccess = public)
        F           % Function (fun) & Proximation (prox) of F
        G           % Function (fun) & Gradient (grad) of G
        H           % Function (fun) & Proximation (prox) of H
        L           % Action (op) & Adjoint (adj) of L
        opt         % Options
    end
    properties (SetAccess = protected)
        x = []
        xOld = []
        y = []
        sigma = 0.0
        tau = 0.0
        theta = 0.0
        res = Inf
        backtrack = []
    end
    methods
        function obj = FBPD(F, G, H, L, opt)
        % Class constructor. 
            if nargin ~= 0
                if nargin < 5; opt = FBPDoptions; end
                obj.F = F;
                obj.G = G;
                obj.H = H;
                obj.L = L;
                obj.opt = opt;
            end            
        end
        
        function obj = updateMaps(obj, F, G, H, L)
            obj.F = F;
            obj.G = G;
            obj.H = H;
            obj.L = L;
        end
        
        function out = minimize(obj, x0)
            obj.initiate(x0);
            k = 0;
            while k < obj.opt.maxIter && obj.res > obj.opt.tol
                obj.iterate;
                obj.backtracking;
                k = k + 1;
            end
            out = obj.x;
        end
    end
    methods (Access = protected)        
        function initiate(obj, x0)
            if obj.opt.resume && ~isempty(obj.y)
                obj.initContinuation(x0);
            else
                obj.fullReset(x0);
            end
        end
        
        function initContinuation(obj, x0)
            obj.x = x0;
            obj.res = Inf;
        end
            
        function fullReset(obj, x0) 
            obj.x = x0;
            obj.y = obj.L.op(x0);
            obj.tau = obj.opt.tau;
            obj.sigma = obj.opt.sigma;
            obj.theta = obj.opt.theta;
            obj.res = Inf;
            obj.backtrack = obj.opt.backtrack;
            obj.initiateBacktrack;
        end
        
        function initiateBacktrack(obj)
            switch class(obj.backtrack)
                case 'BacktrackGLY'
                    obj.backtrack.initParam(obj.tau, obj.sigma, obj.theta);
                    obj.backtrack.initMap(obj.G, obj.L);
                    obj.backtrack.initPoint(obj.x, obj.y);
                case 'BacktrackYH'
                    obj.backtrack.initParam(obj.tau, obj.sigma, obj.theta);
                    obj.backtrack.initMap(obj.G, obj.L);
                    obj.backtrack.initPoint(obj.x, obj.y);
            end
        end
        
        function iterate(obj)
            obj.xOld = obj.x;
            obj.x = obj.x - obj.tau * (obj.G.grad(obj.x) + obj.L.adj(obj.y));
            obj.x = obj.F.prox(obj.x, obj.tau);
            obj.y = obj.y + obj.sigma * obj.L.op((1 + obj.theta) * obj.x ...
                - obj.theta * obj.xOld);
            obj.y = obj.y - obj.sigma *obj.H.prox(obj.y ./ obj.sigma, 1 / obj.sigma);
            obj.res = norm(obj.x(:) - obj.xOld(:)) / norm(obj.xOld(:));
        end
        
        function backtracking(obj)
            if ~isempty(obj.backtrack)
                [obj.tau, obj.sigma] = obj.backtrack.update(obj.x, obj.y);
            end
        end
    end
end
