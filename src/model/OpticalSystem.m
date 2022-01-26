%
%   Parameters of the diffraction model.
%

classdef OpticalSystem
    properties (SetAccess = public)
        m = 0               % number of sampling points
        rM = 0              % distance measurement plane
        mt = 0              % number of rotations
        nk = 0              % grid size in k-space
        grid_factor = 1     % zero padding
        oversampling = 1    % interpolation
        freq = 1            % frequency
        n0 = 1              % refractive index of the medium
        n = [0; 1; 0]       % rotation axis n
        d = 2               % dimension
        lM                  % length of measurement interval/plane
        t0                  % rotation angles
    end
    properties (Dependent)
        nf                  % sampling of f
        nM                  % number of measurement points of u
        xrcv                % grid of the receivers
        k0                  % wave number in surrounding medium
        nu                  % number of measurment points of u per direction
        lk                  % grid size of measurment points of u
        rs                  % maximal size of sample
        rf                  % grid size for f
        r
        z1                  % grid points of f
        z2
        z3
        r1                  % (d-1) dimensional Grid of measurements of u
        r2
        sizeObject          % size(f)
        sizePlane           % size(u)
    end
    methods
        function obj = OpticalSystem(m, mt, nk, rM, d)
        % Class constructor. 
            if nargin ~= 0
                if nargin < 5; d = 2; end
                obj.m = m;
                obj.rM = rM;
                obj.mt = mt;
                obj.nk = nk;
                obj.d = d;
                % Setting default values
                obj.lM = obj.nM / 4; 
%               obj.t0 = (0 : obj.mt - 1) / obj.mt * 2 * pi;
                obj.t0 = mod(-(0:obj.mt-1)/obj.mt*2*pi - pi/2, 2*pi); % angles % pi/2 correction
            end            
        end
        
        function val = get.nf(obj)
            val = obj.grid_factor * obj.oversampling * obj.m; 
        end
        
        function val = get.nM(obj)
            val = obj.nk; 
        end
        
        function val = get.xrcv(obj)
            val = (-obj.nM / 2 : obj.nM / 2 - 1) / obj.nM * 2 * obj.lM;
        end
        
        function val = get.k0(obj)
            val = 2 * pi * obj.freq * obj.n0;
        end
        
        function val = get.nu(obj)
            val = obj.nM;
        end
        
        function val = get.lk(obj)
            val = pi * obj.nM / 2 / obj.lM;
        end
        
        function val = get.rs(obj)
            val = pi * obj.m / (2 * 2 * pi) / sqrt(2);
        end
        
        function val = get.rf(obj)
            val = obj.grid_factor * obj.rs;
        end
        
        function val = get.r(obj)
            val = (-obj.nk / 2 : obj.nk / 2 - 1) / obj.nk * 2 * obj.lk;
        end
        
        function val = get.z1(obj)
            switch obj.d
                case 2
                    [val, ~] = ndgrid((-obj.nf / 2 : obj.nf / 2 - 1) ...
                        / obj.nf * 2 * obj.rf);
                case 3
                    [val, ~, ~] = ndgrid((-obj.nf / 2 : obj.nf / 2 - 1) ...
                        / obj.nf * 2 * obj.rf);
                otherwise
                    error('Invalid dimension');
            end
        end
        
        function val = get.z2(obj)
            switch obj.d
                case 2
                    val = 0*obj.z1;
                case 3
                    [~, val, ~] = ndgrid((-obj.nf / 2 : obj.nf / 2 - 1) ...
                        / obj.nf * 2 * obj.rf);
                otherwise
                    error('Invalid dimension');
            end
        end
        
        function val = get.z3(obj)
            switch obj.d
                case 2
                    [~, val] = ndgrid((-obj.nf / 2 : obj.nf / 2 - 1) ...
                        / obj.nf * 2 * obj.rf);
                case 3
                    [~, ~, val] = ndgrid((-obj.nf / 2 : obj.nf / 2 - 1) ...
                        / obj.nf * 2 * obj.rf);
                otherwise
                    error('Invalid dimension');
            end
        end
        
        function val = get.r1(obj)
            switch obj.d
                case 2
                    [val, ~] = ndgrid(obj.xrcv,0);
                case 3
                    [val, ~] = ndgrid(obj.xrcv,obj.xrcv);
                otherwise
                    error('Invalid dimension');
            end
        end
        
        function val = get.r2(obj)
            switch obj.d
                case 2
                    [~, val] = ndgrid(obj.xrcv,0);
                case 3
                    [~, val] = ndgrid(obj.xrcv,obj.xrcv);
                otherwise
                    error('Invalid dimension');
            end
        end
        
        function val = Green(obj,r)
          switch obj.d
            case 2
              val = 1i/4 * besselh(0,obj.k0*r);
            case 3
              val = exp(1i*obj.k0*r) ./ (4*pi*r);
            otherwise
              error('Invalid dimension');
          end
        end
        
        function out = get.sizeObject(obj)
            switch obj.d
                case 2
                    out = [obj.nf, obj.nf];
                case 3
                    out = [obj.nf, obj.nf, obj.nf];
            end
        end
        
        function out = get.sizePlane(obj)
            switch obj.d
                case 2
                    out = [obj.nu, obj.mt];
                case 3
                    out = [obj.nu, obj.nu, obj.mt];
            end
        end
        
        function obj = set.n(obj,n)
          if isvector(n) && length(n) == 3
            obj.n = n(:);
          elseif  size(n,1) == 3
            obj.n = n;
          else
            error('Rotation axis n must have size 3*1 or 3*mt')
          end
        end
        
        function out = createProbe(obj, type)
            switch lower(type)
                case 'heart1'
                    out = obj.createHeartProbe(obj.m / 6, 1, 1);
                case 'hearts'
                    out = create_model_heart_flo(obj.z1,obj.z3,5,1,1);
                case 'phantom3d1'
                    out = meep_phantom_3d(obj.z1, obj.z2, obj.z3);
                case 'cell3d'
                    out = cell3d(obj.z1, obj.z2, obj.z3, obj.rs);
                case 'ball'
                    out = 1.0 * (obj.z1.^2 + obj.z2.^2 + obj.z3.^2 < .5* obj.rs.^2);
                otherwise
                    error('Unknown probe type');
            end
        end
           
        function out = createHeartProbe(obj, rad, c1, c2)
            [out, ~] = create_model_heart(obj.z1, obj.z3, rad, c1, c2);
        end
    end    
end
