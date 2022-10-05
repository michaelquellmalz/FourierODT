%
%   Forward and backward diffraction.
%

classdef DiffractionModel < handle
    properties (SetAccess = protected)
        os = []
        NFFT = []          
        kappa = []           
        ind = []
        y1 = []
        y2 = []
        y3 = []
        ind_t = []
        FourierGridType = 'uniform'
        FourierGridSize = [];
    end
    properties (SetAccess = public)
        backwardMethod = BackwardCGNE
    end
    methods
        function obj = DiffractionModel(os,gridType)
        % Class constructor. 
            if nargin ~= 0
                if nargin > 1
                  obj.FourierGridType = gridType;
                end
                obj.os = os;
                obj.createLayer;
                obj.setupNFFT;
            end            
        end
        
        function utot = forward(obj, f)
            Ff = obj.objectToLayer(f);
            Ff = Ff .* ((2 * pi)^(-obj.os.d / 2) ...
                * (2 * obj.os.rf / obj.os.nf)^obj.os.d); 
            Fu = zeros(obj.os.sizePlane);
            Fu(obj.ind) = Ff .* sqrt(pi / 2) .* (1i) ...
                .* exp(1i * obj.kappa * obj.os.rM) ./ obj.kappa;
            Fu = Fu ./ ((2 * pi)^(-(obj.os.d-1)/2) ...
                * (2 * obj.os.lM / obj.os.nM)^(obj.os.d-1));
            utot = obj.layerToPlane(Fu) + exp(1i * obj.os.k0 * obj.os.rM);
        end
        
        function utot = forward_convolution(obj,f)
          % only nonzero values of f are important
          z1p = obj.os.z1(f~=0);
          z2p = obj.os.z2(f~=0);
          z3p = obj.os.z3(f~=0);
          fp = f(f~=0);
          uBorn = zeros(obj.os.nM^(obj.os.d-1),obj.os.mt); % initialize
          
          n = obj.os.n .* ones(size(obj.os.t0));
          for t=1:obj.os.mt
            % rotated direction of plane wave uinc
            [s1r,s2r,s3r] = rotation_axis(0,0,1, n(1,t),n(2,t),n(3,t), obj.os.t0(t)); 
            uinc = exp(1i*obj.os.k0*(s1r*z1p+s2r*z2p+s3r*z3p));
            fu = fp .* uinc;
            % rotation of coordinates corresponds to rotation of object f
            [r1r,r2r,r3r] = rotation_axis(...
              obj.os.r1(:),obj.os.r2(:),obj.os.rM,...
              obj.os.n(1,t),obj.os.n(2,t),obj.os.n(3,t),...
              (obj.os.t0(t)));

            % Compute convolution in Matlab
            for j=1:numel(obj.os.r1)
            uBorn(j,t) = sum(obj.os.Green(sqrt((r1r(j)-z1p).^2+(r2r(j)-z2p).^2+(r3r(j)-z3p).^2)) .* fu,'omitnan');
            end
          end

          if (obj.os.d == 3)
            uBorn = reshape(uBorn,obj.os.sizePlane);
          end
          uBorn = uBorn .*(2*obj.os.rf/obj.os.nf)^obj.os.d;
          utot = uBorn + exp(1i * obj.os.k0 * obj.os.rM);
        end
                
        function f = backward(obj, utot, f0)
            if nargin < 3; f0 = zeros(obj.os.sizeObject); end
            uborn = utot - exp(1i * obj.os.k0 * obj.os.rM);
            Fu = obj.planeToLayer(uborn);
            Ff = obj.LayerRemoveWeights(Fu);
            f = obj.LayerToObject(Ff, f0);
        end
                
        function Ff = LayerRemoveWeights(obj, Fu)
            Fu = (2 * pi)^(-(obj.os.d-1)/2)  ...
                * (2 * obj.os.lM / obj.os.nM)^(obj.os.d-1) .* Fu;
            Fu = Fu(obj.ind) .* sqrt(2 / pi) .* (-1i) ...
                .* exp(-1i * obj.kappa * obj.os.rM) .* obj.kappa;
            Ff = Fu(:) ./ ((2 * pi)^(-obj.os.d / 2) ...
                * (2 * obj.os.rf / obj.os.nf)^obj.os.d);
        end
        
        function utot = freeSpaceForward(obj, g)
            u = obj.layerToPlane(g);
            ka = zeros(size(u));
            ka(obj.ind) = obj.kappa;
            u = u ./ exp(-1i * ka * obj.os.rM);
            u = obj.planeToLayer(u);
            utot = u + exp(1i * obj.os.k0 * obj.os.rM);
        end
        
        function g = freeSpaceBackward(obj, utot)
            u = utot - exp(1i * obj.os.k0 * obj.os.rM);
            g = obj.planeToLayer(u);
            ka = zeros(size(g));
            ka(obj.ind) = obj.kappa;
            g = g .* exp(-1i * ka * obj.os.rM);
            g = obj.layerToPlane(g);
        end
        
        function Ff = objectToLayer(obj, f)
            obj.NFFT.fhat = f(:);
            obj.NFFT.nfft_trafo;
            Ff = obj.NFFT.f;
        end
        
        function Fu = planeToLayer(obj, uborn)
            switch obj.os.d
                case 2
                    Fu = fftshift(fft(fftshift(uborn,1), obj.os.nu, 1), 1);
                case 3
                  if 0 % strcmp(obj.FourierGridType, 'uniform')
                    Fu = fftshift(fftshift(fft(fft(fftshift(fftshift(...
                        uborn, 1), 2), obj.os.nu, 1), obj.os.nu, 2), 1), 2);
                  else
                    [k1,k2,~,~] = obj.createFourierGrid3d;
                    k1 = k1(:,:,1);
                    k2 = k2(:,:,1);
                    Fu = obj.planeToLayerNonuniform3d(uborn, k1, k2);
                  end
            end
        end
        
        function Fu = planeToLayerNonuniform3d(obj, uborn, k1, k2)
            Fu = zeros(numel(k1), obj.os.mt);
            p = nfft(2, [obj.os.nu obj.os.nu], prod(obj.FourierGridSize));
            p.x = [k1(:) k2(:)] / (2*obj.os.lk);
            for t = 1:obj.os.mt
              p.fhat = reshape(uborn(:,:,t), [], 1);
              p.nfft_trafo;
              Fu(:,t) = p.f;
            end % for
            Fu = reshape(Fu, size(k1,1), size(k1,2), obj.os.mt);
        end
        
        function w = extractWeights(obj)
            if isempty(obj.backwardMethod.weights)
                obj.computeWeights;
            end
            w = obj.backwardMethod.weights;
        end
    end   
    methods (Access = protected)
        function createLayer(obj)
            switch obj.os.d
                case 2
                    obj.createLayer2d;
                case 3
                    obj.createLayer3d;
            end
        end
        
        function createLayer2d(obj)
            [k1, t] = ndgrid(obj.os.r, obj.os.t0);
            [~, it] = ndgrid(obj.os.r, 1:obj.os.mt);
            obj.FourierGridSize = [obj.os.nk];
            kappa2 = (obj.os.k0^2 - k1.^2);
            obj.ind = kappa2 > 1e-6;
            obj.ind_t = it(obj.ind);
            obj.kappa = sqrt(kappa2(obj.ind));
            y1tmp = k1(obj.ind);
            t = t(obj.ind);
            y3tmp = obj.kappa - obj.os.k0;
            n = reshape(obj.os.n .* ones(size(obj.os.t0)),3,1,obj.os.mt);
            n1 = n(1,1,:) .* ones(obj.os.sizePlane);
            n2 = n(2,1,:) .* ones(obj.os.sizePlane);
            n3 = n(3,1,:) .* ones(obj.os.sizePlane);
            [obj.y1, ~, obj.y3] = rotation_axis(y1tmp, 0, y3tmp, ...
                n1(obj.ind), n2(obj.ind), n3(obj.ind), t);
        end
        
        function createLayer3d(obj)
            obj.FourierGridSize = [obj.os.nk obj.os.nk];
            [k1,k2,t,it] = obj.createFourierGrid3d;
            kappa2 = (obj.os.k0^2 - k1.^2 - k2.^2);
            obj.ind = kappa2>1e-6;
            obj.ind_t = it(obj.ind);
            obj.kappa = sqrt(kappa2(obj.ind));
            y1tmp = k1(obj.ind);
            y2tmp = k2(obj.ind);
            t  =  t(obj.ind);
            y3tmp = obj.kappa - obj.os.k0;
            n = reshape(obj.os.n .* ones(size(obj.os.t0)),3,1,obj.os.mt);
            n1 = n(1,1,:) .* ones(obj.os.sizePlane);
            n2 = n(2,1,:) .* ones(obj.os.sizePlane);
            n3 = n(3,1,:) .* ones(obj.os.sizePlane);
            [obj.y1, obj.y2, obj.y3] = rotation_axis(y1tmp, y2tmp, y3tmp, ...
                n1(obj.ind), n2(obj.ind), n3(obj.ind), t);
        end
        
        function [k1,k2,t,it] = createFourierGrid3d(obj)
          switch obj.FourierGridType
            case 'uniform'
              [k1,k2,t] = ndgrid(obj.os.r, obj.os.r, obj.os.t0);
              [~,~, it] = ndgrid(obj.os.r, obj.os.r, 1:obj.os.mt);
            case 'radial'
              mr = obj.FourierGridSize(1);
              mp = obj.FourierGridSize(2);
              r_line = obj.os.lk*sqrt(2)/(mr+1) * (-mr/2:mr/2-1);
              ph_line = (0:mp-1)*pi/mp;
              [r,ph,t] = ndgrid(r_line, ph_line, obj.os.t0);
              [~,~, it] = ndgrid(r_line, ph_line, 1:obj.os.mt);
              k1 = r.*cos(ph);
              k2 = r.*sin(ph);
            case 'random'
              r1 = rand(obj.FourierGridSize(1),1) * obj.os.lk;
              r2 = rand(obj.FourierGridSize(2),1) * obj.os.lk;
              [k1,k2,t] = ndgrid(r1, r2, obj.os.t0);
              [~,~, it] = ndgrid(r1, r2, 1:obj.os.mt);
            otherwise
              error('Unknown grid type')
          end
        end
        
        function setupNFFT(obj)
            switch obj.os.d
                case 2
                    obj.setupNFFT2d;
                case 3
                    obj.setupNFFT3d;
            end
        end
        
        function setupNFFT2d(obj)
            obj.NFFT = nfft(2, [obj.os.nf, obj.os.nf], numel(obj.y1), ...
                2 * [obj.os.nf, obj.os.nf], 3, ...
                PRE_PHI_HUT + PRE_PSI + NFFT_OMP_BLOCKWISE_ADJOINT, ...
                FFTW_MEASURE);
            obj.NFFT.x = [obj.y1, obj.y3] .* (obj.os.rf / obj.os.nf / pi);
        end
        
        function setupNFFT3d(obj)
            obj.NFFT = nfft(3, [obj.os.nf, obj.os.nf, obj.os.nf], ...
                numel(obj.y1), 1.5 * [obj.os.nf, obj.os.nf, obj.os.nf], 3, ...
                PRE_PHI_HUT + PRE_PSI + NFFT_OMP_BLOCKWISE_ADJOINT,...
                FFTW_MEASURE);
            obj.NFFT.x = [obj.y1, obj.y2, obj.y3] ...
                .* (obj.os.rf / obj.os.nf / pi);
        end
        
        function f = adjointObjectToLayer(obj, Ff)
            obj.NFFT.f = Ff(:);
            obj.NFFT.nfft_adjoint;
            f = reshape(obj.NFFT.fhat, obj.os.nf, obj.os.nf, []);
        end
        
        function uborn = layerToPlane(obj, Fu)
            switch obj.os.d
                case 2
                    uborn = fftshift(ifft(fftshift(Fu, 1), obj.os.nu, 1), 1);
                case 3
                  if 0 %strcmp(obj.FourierGridType, 'uniform')
                    uborn = fftshift(fftshift(ifft(ifft(fftshift(fftshift(...
                        Fu, 1), 2), obj.os.nu, 1), obj.os.nu, 2), 1), 2);
                  else
                    uborn = obj.layerToPlaneNonuniform3d(Fu);
                  end
            end
        end
        
        function uborn = layerToPlaneNonuniform3d(obj, Ff)
          uborn = zeros(obj.os.nu^2, obj.os.mt);
          p = nfft(2, [obj.os.nu obj.os.nu], prod(obj.FourierGridSize));
          [k1,k2,~,~] = obj.createFourierGrid3d;
          k1 = k1(:,:,1);
          k2 = k2(:,:,1);
          p.x = [k1(:) k2(:)] / (2 * obj.os.lk);
          for t = 1:obj.os.mt
            Ff_tmp = Ff(:,:,t);
            p.f = Ff_tmp(:);
            p.nfft_adjoint;
            p.nfft_solver(10);
            uborn(:,t) = p.fhat;
          end % for
            uborn = reshape(uborn,obj.os.nu, obj.os.nu, obj.os.mt);
        end
        
        function f = LayerToObject(obj, Ff, f0)
            switch class(obj.backwardMethod)
                case 'BackwardBP'
                    f = obj.LayerToObjectPB(Ff);
                case 'BackwardCGNE'
                    f = obj.LayerToObjectCGNE(Ff, f0);
                case 'BackwardFBPD'
                    f = obj.LayerToObjectFBPD(Ff, f0);
                otherwise
                    error('Unknown backward method')
            end
        end
        
        function f = LayerToObjectBP(Ff)
            w = extractWeights;
            if obj.backwardMethod.enforceReal
                f = real(obj.adjointObjectToLayer(w .* Ff));
            else
                f = obj.adjointObjectToLayer(w .* Ff);
            end
        end
        
        function f = LayerToObjectCGNE(obj, Ff, f0)
            opt = obj.backwardMethod.CGNEopt;
            w = obj.extractWeights;
            A = obj.CGNE_buildA;
            CG = CGNE(A, Ff, opt);
            f = CG.minimize(f0, w);
        end
        
        function computeWeights(obj)
            switch lower(obj.backwardMethod.weightType)
                case 'uniform'
                    obj.computeWeightsUniform;
                case 'nonuniform'
                    obj.computeWeightsNonuniform;
                otherwise
                    error('Unknown weight type');
            end
        end
        
        function computeWeightsUniform(obj)
            obj.backwardMethod.weights = ones(size(obj.kappa)) / 2 ...
                / (2 * pi)^(obj.os.d) / numel(obj.kappa)...
                * (2 * obj.os.rf / obj.os.nf)^obj.os.d;
        end
        
        function computeWeightsNonuniform(obj)
            switch obj.os.d
                case 2
                    [k1, k2] = ndgrid(obj.os.r, obj.os.t0);
                case 3
                    [k1, k2, ~] = ndgrid(obj.os.r, obj.os.r, obj.os.t0);
            end
            if isempty(obj.os.nd)
            obj.backwardMethod.weights = obj.os.k0 / 2 ...
                * abs(k1(obj.ind)*obj.os.n(2) - k2(obj.ind)*obj.os.n(1)) ...
                ./ (obj.kappa) / (2 * pi)^(obj.os.d) ...
                * obj.volumePlane / numel(obj.kappa)...
                * (2 * obj.os.rf / obj.os.nf)^obj.os.d;
            else
              k1 = k1(obj.ind);
              k2 = k2(obj.ind);
              k3 = obj.kappa - obj.os.k0;
              alpha = repmat(reshape(obj.os.t0,1,1,obj.os.mt), obj.os.nk,obj.os.nk);
              alpha = alpha(obj.ind);
              alphad = ones(size(alpha));
%               n1 = reshape(obj.os.n(1,obj.ind_t), numel(obj.ind_t), 1);
%               n2 = reshape(obj.os.n(2,obj.ind_t), numel(obj.ind_t), 1);
%               n3 = reshape(obj.os.n(3,obj.ind_t), numel(obj.ind_t), 1);
            n = reshape(obj.os.n .* ones(size(obj.os.t0)),3,1,obj.os.mt);
            n1 = n(1,1,:) .* ones(obj.os.sizePlane);
            n2 = n(2,1,:) .* ones(obj.os.sizePlane);
            n3 = n(3,1,:) .* ones(obj.os.sizePlane);
            n1 = n1(obj.ind);
            n2 = n2(obj.ind);
            n3 = n3(obj.ind);
            nd = reshape(obj.os.nd .* ones(size(obj.os.t0)),3,1,obj.os.mt);
            nd1 = nd(1,1,:) .* ones(obj.os.sizePlane);
            nd2 = nd(2,1,:) .* ones(obj.os.sizePlane);
            nd3 = nd(3,1,:) .* ones(obj.os.sizePlane);
            nd1 = nd1(obj.ind);
            nd2 = nd2(obj.ind);
            nd3 = nd3(obj.ind);
%               nd1 = reshape(obj.os.nd(1,obj.ind_t), numel(obj.ind_t), 1);
%               nd2 = reshape(obj.os.nd(2,obj.ind_t), numel(obj.ind_t), 1);
%               nd3 = reshape(obj.os.nd(3,obj.ind_t), numel(obj.ind_t), 1);
              obj.backwardMethod.weights = obj.os.k0 / 2 ...
                ./ (obj.kappa) / (2 * pi)^(obj.os.d) ...
                * obj.volumePlane / numel(obj.kappa)... 
                * (2 * obj.os.rf / obj.os.nf)^obj.os.d ...  
                .* abs( (1-cos(alpha)).*(n3.*(nd1.*k1 ...
                  +nd2.*k2+nd3.*k3)-nd3.*(n1.*k1+n2.*k2+n3.*k3))...
                -n3.*sin(alpha).*(n1.*(nd2.*k3-nd3.*k2))...
                +n2.*(nd3.*k1-nd1.*k3+n3.*(nd1.*k2-nd2.*k1))...
                -alphad.*(n1.*k2-n2.*k1) ...
                +sin(alpha).*(n1.*nd2-n2.*nd1).*(n1.*k1+n2.*k2+n3.*k3));
            end
        end
        
        function val = volumePlane(obj)
          switch obj.os.d
            case 2 
              val = 2 * pi * (2 * obj.os.k0);
            case 3
              val = 2 * pi * (pi * obj.os.k0^2);
          end
        end
        
        function A = CGNE_buildA(obj)
            A.op = @(f) obj.objectToLayer(f);
            if obj.backwardMethod.enforceReal
                A.adj = @(Ff) real(obj.adjointObjectToLayer(Ff));
            else
                A.adj = @(Ff) obj.adjointObjectToLayer(Ff);
            end
        end
        
        function f = LayerToObjectFBPD(obj, Ff, f0)
            opt = obj.backwardMethod.FBPDopt;
            F = obj.FBPD_buildF;
            G = obj.FBPD_buildG(Ff);
            H = obj.FBPD_buildH;
            L = obj.FBPD_buildL;
            PD = obj.setupFBPD(F, G, H, L, opt);
            f = PD.minimize(f0);
        end
        
        function F = FBPD_buildF(obj)
            if obj.backwardMethod.enforceNonneg
                F.prox = @(f, tau) max(0, real(f));
            elseif obj.backwardMethod.enforceReal
                F.prox = @(f, tau) real(f);
            else
                F.prox = @(f, tau) f;
            end
        end
        
        function G = FBPD_buildG(obj, Ff)
            w = obj.extractWeights;
            G.grad = @(f) obj.adjointObjectToLayer(w .* (obj.objectToLayer(f) - Ff));
        end
        
        function H = FBPD_buildH(obj)
            lambda = obj.backwardMethod.lambda;
            H.prox = @(y,gamma) proxL2field(y, gamma*lambda);
        end
        
        function L = FBPD_buildL(~)
            L.op = @(f) gradient(f);
            L.adj = @(y) divergence(y);
        end
        
        function PD = setupFBPD(obj, F, G, H, L, opt)
            if obj.backwardMethod.resume ...
                    && ~isempty(obj.backwardMethod.state)
                PD = obj.backwardMethod.state;
                PD.updateMaps(F, G, H, L);
                PD.opt.resume = true;
            else
                PD = FBPD(F, G, H, L, opt);
                obj.backwardMethod.state = PD;
            end
        end
    end    
end

