classdef CGNE < handle
    properties (SetAccess = public)
        A           % Action (op) & adjoint (adj) of A
        g           % Data
        opt         % Options
    end
    properties (SetAccess = protected)
        x
        xOld
        w
        rk
        pk
        sk
        ak
        bk
        qk
        res
    end
    methods
        function obj = CGNE(A, g, opt)
        % Class constructor. 
            if nargin ~= 0
                if nargin < 3; opt = CGNEoptions; end
                obj.A = A;
                obj.g = g;
                obj.opt = opt;
            end            
        end
        
        function out = minimize(obj, x0, w)
            if nargin < 3; w = ones(size(obj.g)); end
            obj.initiate(x0, w);
            if obj.opt.maxIter == 0
              obj.x = obj.pk; % zero iterations means backpropagation
            end
            k = 0;
            while k < obj.opt.maxIter && obj.res > obj.opt.tol
                obj.iterate;
                k = k + 1;
            end
            out = obj.x;
        end
    end
    methods (Access = protected)
        function initiate(obj, x0, w)
            obj.x = x0;
            obj.xOld = Inf(size(x0));
            obj.w = w;
            obj.rk = obj.g - obj.A.op(x0);
            obj.pk = obj.A.adj(w .* obj.rk);
            obj.sk = obj.pk;
            obj.res = norm(obj.rk, 'fro') / sqrt(length(obj.rk));
        end
        
        function iterate(obj)
            obj.xOld = obj.x;  
            obj.qk = obj.A.op(obj.pk);
            obj.ak = norm(obj.sk(:))^2 ...
                / norm(sqrt(obj.w(:)) .* obj.qk(:))^2;
            obj.x = (obj.x + obj.ak * obj.pk);
            obj.rk = obj.rk - obj.ak * obj.qk;
            nskm1 = norm(obj.sk(:))^2;
            obj.sk = obj.A.adj(obj.w .* obj.rk);
            obj.bk = norm(obj.sk(:))^2 / nskm1;
            obj.pk = obj.sk + obj.bk * obj.pk;
            obj.res = norm(obj.rk(:)) / sqrt(length(obj.rk));
        end
    end
end
