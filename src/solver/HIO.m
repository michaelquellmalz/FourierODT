classdef HIO < handle
    properties (SetAccess = public)
        dm          % Diffraction model
        opt         % Options
    end
    properties (SetAccess = protected)
        f
        fOld
        fInput
        uabs
        res
        supp
    end
    methods
        function obj = HIO(dm, opt)
        % Class constructor. 
            if nargin ~= 0
                if nargin < 2; opt = HIOoptions; end
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
%                  fprintf('It %d,  rk = %g\n', k, obj.res);
            end
            out = obj.enforceObjectConst(obj.f);
        end
    end
    methods (Access = protected)        
        function initiate(obj, utot)
            obj.uabs = abs(utot);
            obj.f = obj.extractF0;
            obj.fInput = obj.f;
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
            u = obj.dm.forward(obj.fInput);
            u = obj.enforcePlaneConst(u);
            obj.f = obj.dm.backward(u, obj.fInput);
            obj.fInput = obj.computeNextInput(obj.f, obj.fOld);
            obj.res = norm(obj.fInput(:) - obj.fOld(:)) / norm(obj.fOld(:));
        end

        function out = enforcePlaneConst(obj, u)
            out = obj.uabs .* exp(1i*angle(u));
        end
        
        function out = computeNextInput(obj, fOut, fIn)          
            switch obj.opt.updateMethod
                case 'HIO'
                    out = obj.HIOupdate(fOut, fIn);
                case 'BIO'
                    out = obj.BIOupdate(fOut, fIn);
                case 'OO'
                    out = obj.OOupdate(fOut);
                otherwise
                    error('Unknown update method');
            end
        end
        
        function out = HIOupdate(obj, fOut, fIn)
            [Df, gamma] = obj.computeUpdate(fOut);
            out = fOut;
            out(gamma) = fIn(gamma) - obj.opt.beta * Df(gamma);
        end
        
        function out = BIOupdate(obj, fOut, fIn)
            [Df, ~] = obj.computeUpdate(fOut);
            out = fIn - obj.opt.beta * Df; 
        end
        
        function out = OOupdate(obj, fOut)
            [Df, ~] = obj.computeUpdate(fOut);
            out = fOut - obj.opt.beta * Df; 
        end
        
        function [Df, gamma] = computeUpdate(obj, fOut)
            [DfRe, gammaRe] = obj.computeRealUpdate(fOut);
            [DfNN, gammaNN] = obj.computeNonnegUpdate(fOut);
            [DfSu, gammaSu] = obj.computeSuppUpdate(fOut);
            Df = DfRe + DfNN;
            Df(gammaSu) = DfSu(gammaSu);
            gamma = gammaRe | gammaNN | gammaSu;
        end
        
        function [Df, gamma] = computeRealUpdate(obj, fOut)
            Df = zeros(size(fOut));
            gamma = false(size(fOut));
            if obj.opt.enforceReal
                gamma = (abs(imag(fOut)) > obj.opt.condTol);
                Df(gamma) = 1i * imag(fOut(gamma));
            end
        end
        
        function [Df, gamma] = computeNonnegUpdate(obj, fOut)
            Df = zeros(size(fOut));
            gamma = false(size(fOut));
            if obj.opt.enforceReal && obj.opt.enforceNonneg
                gamma = (real(fOut) < 0);
                Df(gamma) = real(fOut(gamma));
            end
        end
        
        function [Df, gamma] = computeSuppUpdate(obj, fOut)
            Df = zeros(size(fOut));
            gamma = (abs(fOut) > obj.opt.condTol) & ~obj.supp;
            Df(gamma) = fOut(gamma);
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

