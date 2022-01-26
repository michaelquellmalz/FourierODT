classdef MD < handle
    properties (SetAccess = public)
        dm          % Diffraction model
        opt         % Options
    end
    properties (SetAccess = protected)
        u
        uOld
        uabs
        res
        supp
    end
    methods
        function obj = MD(dm, opt)
        % Class constructor. 
            if nargin ~= 0
                if nargin < 2; opt = DMoptions; end
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
            end
            out = obj.dm.backward(obj.u);
        end
    end
    methods (Access = protected)        
        function initiate(obj, utot)
            obj.uabs = abs(utot);
            obj.u = obj.uabs;
            obj.uOld = Inf(size(obj.u));
            obj.computeSupp;
            obj.res = Inf;
        end
        
        function computeSupp(obj)
          switch obj.dm.os.d
            case 2
              obj.supp = (obj.dm.os.xrcv(:).^2 <= obj.opt.supportRadius^2);
            case 3
              obj.supp = (obj.dm.os.xrcv(:).^2 + obj.dm.os.xrcv(:)'.^2 ...
                  <= obj.opt.supportRadius^2);
          end
        end
        
        function iterate(obj)
            obj.uOld = obj.u;
            g = obj.dm.freeSpaceBackward(obj.u);
            g = obj.enforceObjectConst(g);
            obj.u = obj.dm.freeSpaceForward(g);
            obj.u = obj.enforcePlaneConst(obj.u);
            obj.res = norm(obj.u(:) - obj.uOld(:)) / norm(obj.uOld(:));
        end
        
        function out = enforcePlaneConst(obj, u)
            out = obj.uabs .* exp(1i*angle(u));
        end
        
        function out = enforceObjectConst(obj, g)
            if obj.opt.enforceNonneg
                out = max(0, real(g));
            elseif obj.opt.enforceReal
                out = real(g);
            else
                out = g;
            end
            out = out .* obj.supp;
        end
    end
end

