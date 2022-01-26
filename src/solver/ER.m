classdef ER < handle
    properties (SetAccess = public)
        dm          % Diffraction model
        opt         % Options
    end
    properties (SetAccess = protected)
        f
        fOld
        uabs
        res
        supp
    end
    methods
        function obj = ER(dm, opt)
        % Class constructor. 
            if nargin ~= 0
                if nargin < 2; opt = ERoptions; end
                obj.dm = dm;
                obj.opt = opt;
            end            
        end
        
        function out = traceBack(obj, utot)
            obj.initiate(utot);
            k = 0;
            while k < obj.opt.maxIter && obj.res > obj.opt.tol
                obj.iterate;
                k = k + 1;
                 fprintf('It %d,  rk = %g\n', k, obj.res);
            end
            out = obj.f;
        end
    end
    methods (Access = protected)        
        function initiate(obj, utot)
            obj.uabs = abs(utot);
            obj.f = obj.extractF0;
            obj.fOld = Inf(size(obj.f));
            obj.computeSupp;
            obj.res = Inf;
        end
        
        function computeSupp(obj)
          switch obj.dm.os.d
            case 2
              obj.supp = (obj.dm.os.z1.^2 + obj.dm.os.z3.^2 ...
                  <= obj.opt.supportRadius^2);
            case 3
              obj.supp = (obj.dm.os.z1.^2 + obj.dm.os.z2.^2 + obj.dm.os.z3.^2 ...
                  <= obj.opt.supportRadius^2);
          end
        end
        
        function f = extractF0(obj)
            if isempty(obj.opt.f0)
                f = obj.computeF0;
            else
                f = obj.opt.f0;
            end
        end
        
        function f = computeF0(obj)
            f = obj.dm.backward(obj.uabs);
            f = obj.enforceObjectConst(f);
        end
        
        function iterate(obj)
            obj.fOld = obj.f;
            u = obj.dm.forward(obj.f);
            u = obj.enforcePlaneConst(u);
            obj.f = obj.dm.backward(u, obj.f);
            obj.f = obj.enforceObjectConst(obj.f);
            obj.res = norm(obj.f(:) - obj.fOld(:)) / norm(obj.fOld(:));
        end
        
        function out = enforcePlaneConst(obj, u)
            out = obj.uabs .* exp(1i*angle(u));
        end
        
        function out = enforceObjectConst(obj, f)
            if obj.opt.enforceNonneg
                out = max(0, real(f));
            elseif obj.opt.enforceReal
                out = real(f);
            else
                out = f;
            end
            out(~obj.supp) = 0;
        end
    end
end

